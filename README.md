# GPR Hyperbola Fitting Tool

A human-in-the-loop tool for estimating the dielectric constant (relative permittivity) and depth of buried objects from Ground-Penetrating Radar (GPR) B-scan images.

You annotate the hyperbolic signatures in a scan yourself, and the tool fits each one, works out the wave velocity, and converts that into the permittivity and depth of the object. You stay in control of which signatures get analysed, which matters for noisy or cluttered scans where fully automated detection is not yet reliable.

## How it works

1. Point the tool at a folder of B-scan images.
2. It opens LabelMe so you can trace each hyperbola as a linestrip (from one arm tip, through the apex, to the other arm tip).
3. When you close LabelMe, the tool reads your annotations, fits the GPR hyperbola equation to each one with the Levenberg-Marquardt method, and computes the permittivity and depth.
4. For every image it saves a result figure with the fitted curves drawn on top, plus one CSV summarising every hyperbola.

The permittivity comes from the fitted velocity as epsilon_r = (c / v)^2, and the depth from depth = (t0 / 2) * v.

## Requirements

- Python 3.8 or newer
- numpy
- pillow
- scipy
- matplotlib
- labelme

Install everything with:

```
pip install numpy pillow scipy matplotlib labelme
```

### About LabelMe

The tool uses LabelMe for the annotation step. LabelMe is an open-source image annotation tool written in Python. If the command above did not already install it, you can install it on its own:

```
pip install labelme
```

You do not need to open LabelMe yourself, the script launches it for you. Its documentation and source are here: https://github.com/wkentaro/labelme

## Configuration before you run

Open `gpr_annotate_fit.py` and set the values at the top of the file.

Folder paths:

- `IMAGE_FOLDER`: the folder containing your B-scan images
- `JSON_FOLDER`: where LabelMe will save the annotation JSON files
- `OUTPUT_FOLDER`: where the result figures and CSV are written

Calibration constants. These are specific to your instrument and to how the image was processed, so you need to set them for your own data:

- `DELTA_X_M`: horizontal spacing in metres per pixel (along the survey direction)
- `DELTA_T_NS`: vertical spacing in nanoseconds per pixel (along the time axis)
- `C_M_PER_NS`: speed of light in metres per nanosecond, which you normally leave as it is

The horizontal and vertical spacings come from your acquisition geometry. The values currently in the file are the ones used for the dataset this tool was built on, and they will almost certainly be wrong for a different instrument or a differently processed image, so check them against your own acquisition parameters before trusting the output.

## Running it

```
python gpr_annotate_fit.py
```

LabelMe opens with your image folder loaded. Annotate each hyperbola as a linestrip, save with Ctrl+S, use the arrow buttons to move between images, and close LabelMe once you are done with all of them. The script then processes everything and writes the figures and the CSV to your output folder.

## Output

- One result figure per image, showing the original B-scan with every fitted hyperbola drawn over it, each apex marked, and the permittivity and depth labelled.
- A single `results_summary.csv` listing, for every hyperbola, the apex position, depth, fitted velocity, permittivity, and whether the fit converged.

## Try it yourself

If you want to try the tool on real data, there is a set of GPR B-scan images here:

https://github.com/VatsalTrivedi24/gpr_bscan_image_dataset

Grab a few of those scans, point `IMAGE_FOLDER` at them, set the calibration constants, and start annotating.

## Notes

Each hyperbola needs at least four annotated points for the fit to run. Fits that come out with an unphysical velocity or a permittivity outside the expected range are marked as failed in the CSV instead of being reported, so a failed fit usually means the annotation was too short, too noisy, or did not actually follow a hyperbola.