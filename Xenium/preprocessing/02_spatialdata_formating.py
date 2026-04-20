"""
Convert a Xenium `outs` folder to SpatialData and write as .zarr.

Two modes
---------
Direct (no --split_csv):
    Read the sample folder and write a single .zarr named after --sample_name.

    python 02_spatialdata_formating.py \
        --sample_dir  /path/to/outs \
        --sample_name S21-029273 \
        --result_dir  /path/to/zarr_out

Split (with --split_csv):
    Read a polygon-coordinates CSV (Xenium Explorer selection export) and write
    one .zarr per Selection label found in the CSV.

    python 02_spatialdata_formating.py \
        --sample_dir  /path/to/outs \
        --result_dir  /path/to/zarr_out \
        --split_csv   /path/to/coordinates.csv \
        [--pixel_size 0.2125]
"""

import argparse
import os

import pandas as pd
import spatialdata as sd
import spatialdata_io
from shapely import Polygon


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def region_subsetting(region_df: pd.DataFrame, sdata: sd.SpatialData,
                      pixel_size: float = 0.2125) -> sd.SpatialData:
    """Crop *sdata* to the polygon described by *region_df*.

    Parameters
    ----------
    region_df:
        Rows from the coordinates CSV that belong to a single Selection label.
        Expected columns: 'Selection', 'X' (or first col), 'Y' (or second col).
    sdata:
        Full SpatialData object for this Xenium region.
    pixel_size:
        Xenium pixel size in µm (default 0.2125).
    """
    coords = region_df.set_index('Selection')
    polygon = Polygon(coords.values / pixel_size)
    return sd.polygon_query(sdata, polygon, target_coordinate_system="global")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Format a Xenium outs folder to SpatialData .zarr"
    )
    p.add_argument(
        "--sample_dir", required=True,
        help="Path to the Xenium outs/ directory for this sample/region."
    )
    p.add_argument(
        "--result_dir", required=True,
        help="Output directory where .zarr store(s) will be written."
    )
    p.add_argument(
        "--sample_name", default=None,
        help=(
            "Sample name used as the .zarr filename (direct mode). "
            "Required when --split_csv is NOT provided."
        ),
    )
    p.add_argument(
        "--split_csv", default=None,
        help=(
            "Path to a coordinates CSV exported from Xenium Explorer "
            "(polygon selections). When supplied, one .zarr is written per "
            "Selection label found in the CSV."
        ),
    )
    p.add_argument(
        "--pixel_size", type=float, default=0.2125,
        help="Xenium pixel size in µm used when converting polygon coordinates "
             "(default: 0.2125)."
    )
    return p.parse_args()


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main() -> None:
    args = parse_args()

    os.makedirs(args.result_dir, exist_ok=True)

    print(f"Loading SpatialData from: {args.sample_dir}")
    sdata = spatialdata_io.xenium(args.sample_dir)

    # ------------------------------------------------------------------
    # Split mode: polygon-query each Selection from the coordinates CSV
    # ------------------------------------------------------------------
    if args.split_csv is not None:
        df_split = pd.read_csv(args.split_csv, skiprows=2)
        selections = df_split['Selection'].unique().tolist()
        print(f"Found {len(selections)} selection(s) in {args.split_csv}: {selections}")

        for sel in selections:
            print(f"  Processing selection: {sel}")
            cropped = region_subsetting(
                df_split[df_split['Selection'] == sel],
                sdata,
                pixel_size=args.pixel_size,
            )
            cropped.table.obs['sample'] = sel
            out_path = os.path.join(args.result_dir, f"{sel}.zarr")
            print(f"  Writing → {out_path}")
            cropped.write(out_path)

    # ------------------------------------------------------------------
    # Direct mode: write the full region as a single .zarr
    # ------------------------------------------------------------------
    else:
        if args.sample_name is None:
            raise ValueError(
                "--sample_name is required when --split_csv is not provided."
            )
        sdata.table.obs['sample'] = args.sample_name
        out_path = os.path.join(args.result_dir, f"{args.sample_name}.zarr")
        print(f"Writing → {out_path}")
        sdata.write(out_path)

    print("Done.")


if __name__ == "__main__":
    main()
