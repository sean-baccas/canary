import argparse
import os
import shutil
import subprocess
import sys

# ---------------------------------------------------------------------------
# Geometry, in millimetres, straight from Fig. 1 of the problem statement.
# ---------------------------------------------------------------------------
T = 3.2
PLATE_HX, STEEL_HY, STEEL_HZ = 1.6, 25.0, 63.2
FLANGE_X0, WEB_X0, WEB_X1 = 2.1, 122.1, 125.3
FLANGE_Z0 = 60.0
COIL_A, COIL_RI, COIL_RO, COIL_HZ = 50.0, 25.0, 50.0, 50.0
S2_X = 41.8

# Fig. 2 asks for the average flux density through three sections: S1 (Bz in
# the centre plate at z = 0), S2 (Bx in the top flange at x = 41.8) and S3 (Bz
# in the web at z = 0).  As for TEAM 20, each is reported as the average over a
# thin band of steel carved out of the part.  S1 and S3 lie on the z = 0
# symmetry plane, so their bands run from z = 0 to z = MEASURE_BAND; the S2 band
# is centred on x = S2_X.
MEASURE_BAND = 2.0

MM = 1.0e-3

# Comparison tolerance in millimetres, used for locating entities by position.
# The smallest feature in the model is the 0.5 mm gap between the centre plate
# and the flange tips, so this is still far tighter than anything it must tell
# apart.
TOL = 1.0e-4

def build(args):
    import gmsh

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1 if args.verbose else 0)
    gmsh.model.add("team10")
    occ = gmsh.model.occ

    # Only the eighth of the model with x, y, z >= 0 is meshed.
    air_x, air_y, air_z = args.air_x, args.air_y, args.air_z

    # -- the steel ----------------------------------------------------------
    # Half the centre plate, then one channel cut down to the eighth: the top
    # flange (which also covers the corner) and the web below it.
    plate = occ.addBox(0.0, 0.0, 0.0, PLATE_HX, STEEL_HY, STEEL_HZ)
    flange = occ.addBox(FLANGE_X0, 0.0, FLANGE_Z0, WEB_X1 - FLANGE_X0, STEEL_HY, T)
    web = occ.addBox(WEB_X0, 0.0, 0.0, T, STEEL_HY, FLANGE_Z0)
    # Measurement bands, carved back out of the steel they lie in.  S1 and S3
    # sit on the z = 0 plane; S2 is centred on x = S2_X in the flange.
    s1_band = occ.addBox(0.0, 0.0, 0.0, PLATE_HX, STEEL_HY, MEASURE_BAND)
    s2_band = occ.addBox(S2_X - MEASURE_BAND / 2, 0.0, FLANGE_Z0, MEASURE_BAND, STEEL_HY, T)
    s3_band = occ.addBox(WEB_X0, 0.0, 0.0, T, STEEL_HY, MEASURE_BAND)

    # -- the coil -----------------------------------------------------------
    # A rounded square ring: a filleted box with a smaller filleted box cut out
    # of it.  Filleting the four vertical edges of a box of half-width a with
    # radius r gives a rounded square whose corner arcs are centred at
    # (+-(a - r), +-(a - r)), so both boxes share the arc centres (+-50, +-50)
    # that Fig. 1 draws.  The ring is then cut down to the eighth being meshed.
    def rounded_box(half, radius, z0, dz):
        tag = occ.addBox(-half, -half, z0, 2 * half, 2 * half, dz)
        occ.synchronize()
        curves = set()
        for face in gmsh.model.getBoundary([(3, tag)], combined=False, oriented=False):
            for c in gmsh.model.getBoundary([face], combined=False, oriented=False):
                curves.add(abs(c[1]))
        vertical = []
        for c in curves:
            x0, y0, zc0, x1, y1, zc1 = gmsh.model.getBoundingBox(1, c)
            # OCC pads bounding boxes by its own gap tolerance, so this test is
            # loose compared with 1e-9 but still far tighter than any feature.
            if max(x1 - x0, y1 - y0) < TOL and abs(zc1 - zc0 - dz) < TOL:
                vertical.append(c)
        return occ.fillet([tag], vertical, [radius])[0][1]

    outer = rounded_box(COIL_A + COIL_RO, COIL_RO, 0.0, COIL_HZ)
    inner = rounded_box(COIL_A + COIL_RI, COIL_RI, -5.0, COIL_HZ + 10.0)
    ring = occ.cut([(3, outer)], [(3, inner)])[0][0][1]
    eighth = occ.addBox(0.0, 0.0, -5.0, 200.0, 200.0, COIL_HZ + 10.0)
    coil = occ.intersect([(3, ring)], [(3, eighth)])[0][0][1]

    # -- free space and the conformal assembly ------------------------------
    box = occ.addBox(0.0, 0.0, 0.0, air_x, air_y, air_z)
    solids = [plate, flange, web, s1_band, s2_band, s3_band, coil]
    _, fmap = occ.fragment([(3, box)], [(3, t) for t in solids])
    occ.synchronize()

    # -- name the volumes ---------------------------------------------------
    # fragment reports, for each input solid, the pieces it ended up as.  The
    # first entry is the free space box, which holds every piece of the model;
    # whatever is not claimed by a solid is the air.
    pieces = [set(t for _, t in ents) for ents in fmap]
    plate_tags = pieces[1]
    channel_tags = pieces[2] | pieces[3]
    s1_tags, s2_tags, s3_tags, coil_tags = pieces[4], pieces[5], pieces[6], pieces[7]
    air_tags = pieces[0] - plate_tags - channel_tags - coil_tags
    # The measurement bands overlap the plate and the channel, so take them back out.
    plate_tags -= s1_tags
    channel_tags -= s2_tags | s3_tags

    named = [
        ("Plate", plate_tags),
        ("PlateS1", s1_tags),
        ("Channel", channel_tags),
        ("ChannelS2", s2_tags),
        ("ChannelS3", s3_tags),
        ("Coil", coil_tags),
        ("Air", air_tags),
    ]
    for name, tags in named:
        if not tags:
            sys.exit(f"no volumes were assigned to block {name}")
        gmsh.model.setPhysicalName(3, gmsh.model.addPhysicalGroup(3, sorted(tags)), name)
    print("volumes: " + ", ".join(f"{n}={len(t)}" for n, t in named))

    # -- name the surfaces --------------------------------------------------
    # x = 0 and y = 0 contain the coil axis and carry n x A = 0, as does the far
    # field.  z = 0 is the mid-plane of the coil, where B is normal to the plane;
    # it carries the natural condition, so it must be kept apart from the others.
    surf = {}

    def add(name, tag):
        surf.setdefault(name, []).append(tag)

    for _, tag in gmsh.model.getEntities(2):
        x0, y0, z0, x1, y1, z1 = gmsh.model.getBoundingBox(2, tag)
        if abs(x1 - x0) < TOL and abs(x0) < TOL:
            add("symmetry_x", tag)
        elif abs(y1 - y0) < TOL and abs(y0) < TOL:
            add("symmetry_y", tag)
        elif abs(z1 - z0) < TOL and abs(z0) < TOL:
            add("symmetry_z", tag)
        elif (abs(x1 - x0) < TOL and abs(x0 - air_x) < TOL) or \
             (abs(y1 - y0) < TOL and abs(y0 - air_y) < TOL) or \
             (abs(z1 - z0) < TOL and abs(z0 - air_z) < TOL):
            add("far_field", tag)

    for name, tags in sorted(surf.items()):
        gmsh.model.setPhysicalName(2, gmsh.model.addPhysicalGroup(2, tags), name)
    print("surfaces: " + ", ".join(f"{k}={len(v)}" for k, v in sorted(surf.items())))
    missing = [n for n in ("far_field", "symmetry_x", "symmetry_y", "symmetry_z")
               if n not in surf]
    if missing:
        sys.exit("missing expected surface groups: " + ", ".join(missing))

    # -- element sizing -----------------------------------------------------
    # Box fields, one per feature.  Each grades from its own size inside the box
    # to the far-field size over the transition thickness outside it.  The steel
    # is given one box per part rather than one round the whole channel, which
    # would also fill the air inside the channel window with steel-sized elements.
    def box_field(bounds, size, thickness):
        f = gmsh.model.mesh.field.add("Box")
        (x0, y0, z0), (x1, y1, z1) = bounds
        for k, v in (
            ("XMin", x0), ("XMax", x1), ("YMin", y0), ("YMax", y1),
            ("ZMin", z0), ("ZMax", z1), ("VIn", size), ("VOut", args.air_size),
            ("Thickness", thickness),
        ):
            gmsh.model.mesh.field.setNumber(f, k, v)
        return f

    # Isotropic tets, so the element size in the steel is what sets the number
    # of elements through its thickness.
    steel_size = T / args.through_thickness
    fields = [
        # the 0.5 mm gap between the centre plate and the flange tip
        box_field(((PLATE_HX - 1.0, 0.0, FLANGE_Z0 - 1.0),
                   (FLANGE_X0 + 1.0, STEEL_HY + 1.0, STEEL_HZ + 1.0)),
                  args.gap_size, 3.0),
        # the steel: plate, flange and web
        box_field(((0.0, 0.0, 0.0), (PLATE_HX, STEEL_HY, STEEL_HZ)), steel_size, 5.0),
        box_field(((FLANGE_X0, 0.0, FLANGE_Z0), (WEB_X1, STEEL_HY, STEEL_HZ)), steel_size, 5.0),
        box_field(((WEB_X0, 0.0, 0.0), (WEB_X1, STEEL_HY, FLANGE_Z0)), steel_size, 5.0),
        # the coil and the air immediately around the device
        box_field(((0.0, 0.0, 0.0),
                   (WEB_X1 + 5.0, COIL_A + COIL_RO + 5.0, STEEL_HZ + 5.0)),
                  args.coil_size, 40.0),
    ]
    fmin = gmsh.model.mesh.field.add("Min")
    gmsh.model.mesh.field.setNumbers(fmin, "FieldsList", fields)
    gmsh.model.mesh.field.setAsBackgroundMesh(fmin)
    for opt, val in (
        ("Mesh.MeshSizeExtendFromBoundary", 0),
        ("Mesh.MeshSizeFromPoints", 0),
        ("Mesh.MeshSizeFromCurvature", 0),
        ("Mesh.MeshSizeMax", args.air_size),
        ("Mesh.Algorithm3D", 10),  # HXT
        ("Mesh.Optimize", 1),
        ("Mesh.OptimizeNetgen", 0),
    ):
        gmsh.option.setNumber(opt, val)

    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.ScalingFactor", MM)  # millimetres -> metres
    gmsh.write(args.msh)
    n_tet = len(gmsh.model.mesh.getElementsByType(4)[0])
    gmsh.finalize()
    print(f"wrote {args.msh}: {n_tet} tetrahedra")
    return n_tet

def to_exodus(args):
    """Convert to Exodus so the block and sideset names reach MFEM.

    MFEM picks up named attribute sets from Exodus only; its gmsh reader gives
    bare integer attributes.  MOOSE reads the physical names out of the .msh and
    writes them as element block and sideset names, so a --mesh-only pass is all
    that is needed.
    """
    exe = args.exe or shutil.which("moose_test-opt")
    if not exe:
        moose_dir = os.environ.get("MOOSE_DIR")
        if moose_dir:
            cand = os.path.join(moose_dir, "test", "moose_test-opt")
            exe = cand if os.path.exists(cand) else None
    if not exe:
        sys.exit("no MOOSE executable found; pass --exe or set MOOSE_DIR")

    here = os.path.dirname(os.path.abspath(args.msh)) or "."
    conv = os.path.join(here, ".team20_convert.i")
    with open(conv, "w") as fh:
        fh.write(
            "# Written by team10_mesh.py; converts the gmsh mesh to Exodus so that the\n"
            "# block and sideset names reach MFEM, which reads named attribute sets from\n"
            "# Exodus only.\n"
            "[Mesh]\n  [fmg]\n    type = FileMeshGenerator\n"
            f"    file = {os.path.basename(args.msh)}\n  []\n[]\n"
            "[Problem]\n  solve = false\n[]\n[Executioner]\n  type = Steady\n[]\n"
        )
    try:
        subprocess.run(
            [exe, "-i", os.path.basename(conv), "--mesh-only", os.path.basename(args.exodus)],
            cwd=here, check=True,
            stdout=None if args.verbose else subprocess.DEVNULL,
        )
    finally:
        os.path.exists(conv) and os.remove(conv)
    print(f"wrote {args.exodus}")

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gap-size", type=float, default=0.25,
                   help="element size in the 0.5 mm gap between the centre plate "
                        "and the flange tips (mm)")
    p.add_argument("--through-thickness", type=int, default=4,
                   help="number of elements through the 3.2 mm steel; sets the "
                        "steel element size to 3.2 / N")
    p.add_argument("--coil-size", type=float, default=6.0,
                   help="element size in the coil and the air around the device (mm)")
    p.add_argument("--air-size", type=float, default=40.0,
                   help="element size at the outer boundary (mm)")
    p.add_argument("--air-x", type=float, default=350.0,
                   help="extent of the free space box in x (mm)")
    p.add_argument("--air-y", type=float, default=300.0,
                   help="extent of the free space box in y (mm)")
    p.add_argument("--air-z", type=float, default=300.0,
                   help="extent of the free space box in z (mm)")
    p.add_argument("--msh", default="team10.msh")
    p.add_argument("--exodus", default="team10.e")
    p.add_argument("--exe", help="MOOSE executable used for the Exodus conversion")
    p.add_argument("--no-exodus", action="store_true", help="stop after writing the .msh")
    p.add_argument("-v", "--verbose", action="store_true")
    args = p.parse_args()

    if args.through_thickness < 1:
        p.error("--through-thickness must be at least 1")

    build(args)
    if not args.no_exodus:
        to_exodus(args)


if __name__ == "__main__":
    main()
