import ROOT
from ROOT import TGraph, TCanvas, TMath, TF1, TLegend

from showerrec import Seg  # Import the ShowerData class

def HuberWeight(residual, robustness):
    threshold = robustness
    if abs(residual) <= threshold:
        return 1.0
    else:
        return threshold / abs(residual)

def robust_axis(axis, segments, n_segments, robustness):
    max_iterations = 10
    convergence_threshold = 1e-3

    x_mean, y_mean, z_mean = 0.0, 0.0, 0.0
    tx_mean, ty_mean = 0.0, 0.0

    for i in range(n_segments):
        segment = segments.At(i)
        x_mean += segment.eX
        y_mean += segment.eY
        z_mean += segment.eZ
        tx_mean += segment.eTX
        ty_mean += segment.eTY

    x_mean /= n_segments
    y_mean /= n_segments
    z_mean /= n_segments
    tx_mean /= n_segments
    ty_mean /= n_segments

    for _ in range(max_iterations):
        x_sum, y_sum, z_sum = 0.0, 0.0, 0.0
        tx_sum, ty_sum = 0.0, 0.0
        weight_sum = 0.0

        for i in range(n_segments):
            segment = segments.At(i)

            dx = segment.eX - x_mean
            dy = segment.eY - y_mean
            dz = segment.eZ - z_mean
            dtx = segment.eTX - tx_mean
            dty = segment.eTY - ty_mean

            weight_x = HuberWeight(dx, robustness)
            weight_y = HuberWeight(dy, robustness)
            weight_z = HuberWeight(dz, robustness)
            weight_tx = HuberWeight(dtx, robustness)
            weight_ty = HuberWeight(dty, robustness)

            weight = (weight_x + weight_y + weight_z + weight_tx + weight_ty) / 5.0

            x_sum += segment.eX * weight
            y_sum += segment.eY * weight
            z_sum += segment.eZ * weight
            tx_sum += segment.eTX * weight
            ty_sum += segment.eTY * weight
            weight_sum += weight

        new_x_mean = x_sum / weight_sum
        new_y_mean = y_sum / weight_sum
        new_z_mean = z_sum / weight_sum
        new_tx_mean = tx_sum / weight_sum
        new_ty_mean = ty_sum / weight_sum

        if (abs(new_x_mean - x_mean) < convergence_threshold and
            abs(new_y_mean - y_mean) < convergence_threshold and
            abs(new_z_mean - z_mean) < convergence_threshold and
            abs(new_tx_mean - tx_mean) < convergence_threshold and
            abs(new_ty_mean - ty_mean) < convergence_threshold):
            break

        x_mean, y_mean, z_mean = new_x_mean, new_y_mean, new_z_mean
        tx_mean, ty_mean = new_tx_mean, new_ty_mean

    axis.set_attributes(x_mean, y_mean, z_mean, tx_mean, ty_mean, axis.pid)


def plot_robust_profile(axis_collections, segments, canvas=None, filename="shower.pdf"):
    """
    Plots the robust fit profile of the shower axis.
    
    Parameters:
    - axis_collections: List of axis parameters. Each entry is a list or tuple [X0, Y0, Z0, TX, TY].
    - segments: A ROOT TObjArray containing the segments of the shower.
    - canvas: A ROOT TCanvas for plotting. If None, a new canvas will be created.
    - filename: Name of the output PDF file to save the plots.
    """
    # Create a canvas if none is provided
    if canvas is None:
        canvas = TCanvas("canvas", "Robust Fit Profile", 800, 600)
    
    canvas.Clear()
    canvas.Divide(2, 2)  # Divide canvas into 2x2 sub-pads

    # Initialize graphs for X-Z and Y-Z projections
    grx, gry = TGraph(), TGraph()
    for i in range(segments.GetEntriesFast()):
        segment = segments.At(i)
        grx.SetPoint(grx.GetN(), segment.eZ, segment.eX)
        gry.SetPoint(gry.GetN(), segment.eZ, segment.eY)

    # Determine the Z range for the plots
    zmin = TMath.MinElement(grx.GetN(), grx.GetX())
    zmax = TMath.MaxElement(grx.GetN(), grx.GetX())

    # Draw the X-Z projection
    canvas.cd(1)
    grx.Draw("AP")  # "A" for axis, "P" for points
    grx.SetTitle("X-Z Projection")
    grx.GetXaxis().SetTitle("Z")
    grx.GetYaxis().SetTitle("X")
    legx = TLegend(0.1, 0.7, 0.48, 0.9)  # Legend for X-Z projection

    # Draw the Y-Z projection
    canvas.cd(2)
    gry.Draw("AP")
    gry.SetTitle("Y-Z Projection")
    gry.GetXaxis().SetTitle("Z")
    gry.GetYaxis().SetTitle("Y")
    legy = TLegend(0.1, 0.7, 0.48, 0.9)  # Legend for Y-Z projection

    # Plot the axis fits for each set of parameters in axis_collections
    for i, axis in enumerate(axis_collections):
        x0, y0, z0, tx, ty = axis

        # Define and configure the X-Z fit function
        fit_xz = TF1(f"fitXZ_{i}", "[0] + [1] * (x - [2])", zmin, zmax)  # X = X0 + TX * (Z - Z0)
        fit_xz.SetParameters(x0, tx, z0)
        fit_xz.SetLineColor(ROOT.kBlue if i == 0 else ROOT.kRed + i)
        fit_xz.SetLineWidth(2)

        # Draw the X-Z fit on the graph
        canvas.cd(1)
        fit_xz.Draw("L SAME")  # "L" for line
        label = "1-segment axis" if i == 0 else f"ROB={0.01 * (i + 1):.2f}"
        legx.AddEntry(fit_xz, label, "l")  # Add entry to the legend

        # Define and configure the Y-Z fit function
        fit_yz = TF1(f"fitYZ_{i}", "[0] + [1] * (x - [2])", zmin, zmax)  # Y = Y0 + TY * (Z - Z0)
        fit_yz.SetParameters(y0, ty, z0)
        fit_yz.SetLineColor(ROOT.kBlue if i == 0 else ROOT.kRed + i)
        fit_yz.SetLineWidth(2)

        # Draw the Y-Z fit on the graph
        canvas.cd(2)
        fit_yz.Draw("L SAME")
        legy.AddEntry(fit_yz, label, "l")  # Add entry to the legend

    # Finalize the X-Z projection legend
    canvas.cd(1)
    legx.SetHeader("X-Z Fit Axes", "C")  # Centered header
    legx.Draw()

    # Finalize the Y-Z projection legend
    canvas.cd(2)
    legy.SetHeader("Y-Z Fit Axes", "C")  # Centered header
    legy.Draw()

    # Save the canvas to a file
    canvas.Print(filename)
