import ROOT
from ROOT import TGraph, TCanvas, TMath, TF1, TLegend, TLine, TH1I, TH1D, TText, TNtuple
import pandas as pd
import numpy as np
from collections import defaultdict

class Seg:
    pid: int
    def __init__(self):
        self.X = 0.0
        self.Y = 0.0
        self.Z = 0.0
        self.TX = 0.0
        self.TY = 0.0
        self.E0 = 0.0
        self.pid = 0

    def set_attributes(self, X=0.0, Y=0.0, Z=0.0, TX=0.0, TY=0.0, pid=0):
        self.X = X
        self.Y = Y
        self.Z = Z
        self.TX = TX
        self.TY = TY
        self.pid = pid

    def propagate_to(self, z):
        self.X = self.X + self.TX * (z-self.Z)
        self.Y = self.Y + self.TY * (z-self.Z)
        self.Z = z

def calc_dmin(seg1, seg2, dminz=None):
    if seg1.pid == seg2.pid:
        return 0.0

    x1, y1, z1, ax1, ay1 = seg1.X, seg1.Y, seg1.Z, seg1.TX, seg1.TY
    x2, y2, z2, ax2, ay2 = seg2.X, seg2.Y, seg2.Z, seg2.TX, seg2.TY

    s1bunsi = (ax2 * ax2 + ay2 * ay2 + 1) * (ax1 * (x2 - x1) + ay1 * (y2 - y1) + z2 - z1) - \
              (ax1 * ax2 + ay1 * ay2 + 1) * (ax2 * (x2 - x1) + ay2 * (y2 - y1) + z2 - z1)
    s1bunbo = (ax1 * ax1 + ay1 * ay1 + 1) * (ax2 * ax2 + ay2 * ay2 + 1) - \
              (ax1 * ax2 + ay1 * ay2 + 1) * (ax1 * ax2 + ay1 * ay2 + 1)
    s2bunsi = (ax1 * ax2 + ay1 * ay2 + 1) * (ax1 * (x2 - x1) + ay1 * (y2 - y1) + z2 - z1) - \
              (ax1 * ax1 + ay1 * ay1 + 1) * (ax2 * (x2 - x1) + ay2 * (y2 - y1) + z2 - z1)
    s2bunbo = (ax1 * ax1 + ay1 * ay1 + 1) * (ax2 * ax2 + ay2 * ay2 + 1) - \
              (ax1 * ax2 + ay1 * ay2 + 1) * (ax1 * ax2 + ay1 * ay2 + 1)

    s1 = s1bunsi / s1bunbo
    s2 = s2bunsi / s2bunbo

    p1x = x1 + s1 * ax1
    p1y = y1 + s1 * ay1
    p1z = z1 + s1 * 1
    p2x = x2 + s2 * ax2
    p2y = y2 + s2 * ay2
    p2z = z2 + s2 * 1

    p1p2 = np.sqrt((p1x - p2x) ** 2 + (p1y - p2y) ** 2 + (p1z - p2z) ** 2)

    if dminz is not None:
        dminz[0] = z1 - p1z

    return p1p2

def ErecShowerMaximum(showermax):
    return TMath.Power(10, 0.060 * showermax + 1.39)

def ErecNsegSumShowerMax(nsegSumShowerMax):
    p0, p1, p2, p3 = 22.9916, 2.4266, 0.0189621, 8.45026e-5
    x = nsegSumShowerMax
    return p0 + p1 * x + p2 * x**2 + p3 * x**3

def fit_track(t):
    if t.N() <= 1:
        return

    grx, gry = TGraph(), TGraph()

    for i in range(min(t.N(), 4)):
        s = t.GetSegment(i)
        grx.SetPoint(grx.GetN(), s.Z, s.X)
        gry.SetPoint(grx.GetN(), s.Z, s.Y)

    grx.Fit("pol1", "ROB=0.5 R0Q")
    gry.Fit("pol1", "ROB=0.5 R0Q")

    t.SetTX(grx.GetFunction("pol1").GetParameter(1))
    t.SetTY(gry.GetFunction("pol1").GetParameter(1))
    t.PrintNice()

def fit_shower_axis(axis, segments, canvas=None, filename="shower.pdf"):
    canvas.Clear()
    canvas.Divide(2, 2)

    grx, gry = TGraph(), TGraph()
    for i in range(segments.GetEntriesFast()):
        s = segments.At(i)
        grx.SetPoint(grx.GetN(), s.Z(), s.X())
        gry.SetPoint(grx.GetN(), s.Z(), s.Y())

    xmin, xmax = TMath.MinElement(grx.GetN(), grx.GetX()), TMath.MaxElement(grx.GetN(), grx.GetX())
    ymin, ymax = TMath.MinElement(gry.GetN(), gry.GetX()), TMath.MaxElement(gry.GetN(), gry.GetX())

    fits = [("fx", xmin, xmax), ("fy", ymin, ymax)]
    for name, min_val, max_val in fits:
        func = TF1(name, "[0]*x + [1]", min_val, max_val)
        func.SetLineColor(ROOT.kBlue)
        grx.Fit(func, "R0Q")

    canvas.cd(1)
    grx.Draw("AP")
    canvas.cd(2)
    gry.Draw("AP")
    canvas.Print(filename)

def reset_axis(axis, segments, plate):
    grx, gry = TGraph(), TGraph()
    for i in range(plate):
        s = segments.At(i)
        grx.SetPoint(grx.GetN(), s.Z, s.X)
        gry.SetPoint(grx.GetN(), s.Z, s.Y)

    grx.Fit("pol1", "Q")
    gry.Fit("pol1", "Q")

    x = grx.GetFunction("pol1").GetParameter(0)
    y = gry.GetFunction("pol1").GetParameter(0)
    z = 0.0
    tx = grx.GetFunction("pol1").GetParameter(1)
    ty = gry.GetFunction("pol1").GetParameter(1)
    axis.set_attributes(x, y, z, tx, ty)

def initial_axis(segments):
    # use the segments of the primary track on the first 4 plates as an initial axis estimation
    # segments = df["trkID==one_primary_track"]
    max_pid = segments["PID"].max()
    plate = max(4, max_pid)

    segments_new = segments[(segments["PID"]<plate)]
    
    grx, gry = TGraph(), TGraph()
    for _, seg in segments_new.iterrows():
        grx.SetPoint(grx.GetN(), seg["Z"], seg["X"])
        gry.SetPoint(gry.GetN(), seg["Z"], seg["Y"])

    grx.Fit("pol1", "Q")
    gry.Fit("pol1", "Q")

    x = grx.GetFunction("pol1").GetParameter(0)
    y = gry.GetFunction("pol1").GetParameter(0)
    z = 0.0
    tx = grx.GetFunction("pol1").GetParameter(1)
    ty = gry.GetFunction("pol1").GetParameter(1)
    axis = Seg()
    axis.set_attributes(x, y, z, tx, ty)
    return axis


def count_cylinder(axis, segments, z0, d2_max=100**2, dt2_max=0.01**2, dmin_max=50):
    # segments = df

    selected_segments = []
    axis.propagate_to(z0)

    segments_new = segments[(abs(segments["X"]-axis.X)<1000) & (abs(segments["Y"]-axis.Y)<1000) & (segments["Z"]>z0) & (segments["Z"]<z0+80000)]

    for _, seg in segments_new.iterrows():
        axis.propagate_to(seg["Z"])
        x = seg["X"]
        y = seg["Y"]
        # print(f"axis X={axis.X:.1f}, Y={axis.Y:.1f}, Z={axis.Z:.1f}; seg X={x:.1f}, Y={y:.1f}")
        d2 = (seg["X"] - axis.X)**2 + (seg["Y"] - axis.Y)**2
        # print(f"d2 = {d2}")
        if d2 > d2_max:
            continue

        dt2 = (seg["TX"] - axis.TX)**2 + (seg["TY"] - axis.TY)**2
        # print(f"dt2 = {dt2}")
        if dt2 > dt2_max:
            continue

        new_seg = Seg()
        new_seg.set_attributes(seg["X"], seg["Y"], seg["Z"], seg["TX"], seg["TY"], seg["PID"])
        dmin = calc_dmin(axis, new_seg)
        if dmin > dmin_max:
            continue

        selected_segments.append(seg)

    axis.propagate_to(z0)
    return pd.DataFrame(selected_segments) # return the segments in the cylinder

def plot_profile(new_collections, segments, new_segments, nt_profile, ntuple, labels, top_indices, c1=None, t="shower.pdf", best=0):
    if c1 is None:
        c1 = ROOT.TCanvas("c1", "Canvas", 800, 600)
    c1.Clear()
    c1.Divide(2, 2)

    # Plot 1: All segments in one event
    c1.cd(1)
    iplmin0 = int(ntuple.GetMinimum("plate"))
    iplmax0 = int(ntuple.GetMaximum("plate"))
    hProfile0 = TH1D("hProfile0", "All segments in one event", (iplmax0 - iplmin0) + 1, iplmin0 - 0.5, iplmax0 + 0.5)
    hProfile0.SetXTitle("PID (Pattern ID)")
    hProfile0.SetYTitle("Number of segments")
    ntuple.Draw("plate >> hProfile0", "", "goff")
    hProfile0.SetMarkerStyle(20)
    hProfile0.SetMarkerSize(0.5)
    hProfile0.SetLineColor(ROOT.kBlack)

    # Blur profile to reduce statistical fluctuation
    hBlurred0 = hProfile0.Clone("hBlurred0")
    for ibin in range(1, hProfile0.GetNbinsX() + 1):
        sum_val = sum(
            hProfile0.GetBinContent(i) * TMath.Gaus(i, ibin, 3.5, True)
            for i in range(1, hProfile0.GetNbinsX() + 1)
        )
        hBlurred0.SetBinContent(ibin, sum_val)
        hBlurred0.SetBinError(ibin, 0)
        hProfile0.SetBinError(ibin, 0)
    hBlurred0.SetLineColor(ROOT.kBlue)

    hProfile0.Draw("p")
    hBlurred0.Draw("l same")
    maxbin0 = hBlurred0.GetMaximumBin()
    maxval0 = hBlurred0.GetBinContent(maxbin0)
    maxpid0 = hBlurred0.GetBinCenter(maxbin0)

    l0 = TLine()
    l0.SetLineColor(ROOT.kBlue)
    l0.SetLineWidth(3)
    l0.DrawLine(maxpid0, 0, maxpid0, maxval0)
    tx = TText()
    showermaximum0 = maxpid0 - iplmin0 + 1
    tx.DrawText(maxpid0, maxval0 * 0.3, f"shower max at {int(showermaximum0)}")
    print(f"maxbin = {maxbin0}, maxval = {maxval0:.1f}, maxpid = {maxpid0:.1f}, deltaPID = {showermaximum0:.1f}")

    c1.Modified()
    c1.Update()

    # Plot 2: Best-cylindered segments in one event
    c1.cd(2)
    iplmin = int(nt_profile.GetMinimum("plate"))
    iplmax = int(nt_profile.GetMaximum("plate"))
    hProfile = TH1D("hProfile", "Best-cylindered segs in one event", (iplmax - iplmin) + 1, iplmin - 0.5, iplmax + 0.5)
    hProfile.SetXTitle("PID (Pattern ID)")
    hProfile.SetYTitle("Number of segments")
    nt_profile.Draw("plate >> hProfile", "", "goff")
    hProfile.SetMarkerStyle(20)
    hProfile.SetMarkerSize(0.5)
    hProfile.SetLineColor(ROOT.kBlack)

    # Blur profile
    hBlurred = hProfile.Clone("hBlurred")
    for ibin in range(1, hProfile.GetNbinsX() + 1):
        sum_val = sum(
            hProfile.GetBinContent(i) * TMath.Gaus(i, ibin, 3.5, True)
            for i in range(1, hProfile.GetNbinsX() + 1)
        )
        hBlurred.SetBinContent(ibin, sum_val)
        hBlurred.SetBinError(ibin, 0)
        hProfile.SetBinError(ibin, 0)
    hBlurred.SetLineColor(ROOT.kBlue)

    hProfile.Draw("p")
    hBlurred.Draw("l same")
    maxbin = hBlurred.GetMaximumBin()
    maxval = hBlurred.GetBinContent(maxbin)
    maxpid = hBlurred.GetBinCenter(maxbin)

    l = TLine()
    l.SetLineColor(ROOT.kBlue)
    l.SetLineWidth(3)
    l.DrawLine(maxpid, 0, maxpid, maxval)
    showermaximum = maxpid - iplmin + 1
    tx.DrawText(maxpid, maxval * 0.4, f"shower max at {int(showermaximum)}")
    tx.DrawText(maxpid, maxval * 0.3, f"with axis fitted by {labels[best]}")
    print(f"maxbin = {maxbin}, maxval = {maxval:.1f}, maxpid = {maxpid:.1f}, deltaPID = {showermaximum:.1f}")

    c1.Modified()
    c1.Update()

    # Plot 3 and 4: Segments and newSegments
    grx0, gry0 = TGraph(), TGraph()
    for i, seg in enumerate(segments):
        grx0.SetPoint(grx0.GetN(), seg.eZ, seg.eX)
        gry0.SetPoint(gry0.GetN(), seg.eZ, seg.eY)

    grx, gry = TGraph(), TGraph()
    for i, seg in enumerate(new_segments):
        grx.SetPoint(grx.GetN(), seg.eZ, seg.eX)
        gry.SetPoint(gry.GetN(), seg.eZ, seg.eY)

    zmin = min(grx0.GetX())
    zmax = max(grx0.GetX())

    c1.cd(3)
    grx0.SetMarkerStyle(20)
    grx0.SetMarkerSize(0.5)
    grx0.GetXaxis().SetTitle("Z/um")
    grx0.GetYaxis().SetTitle("X/um")
    grx0.Draw("ap")
    grx.SetMarkerStyle(20)
    grx.SetMarkerSize(0.6)
    grx.SetMarkerColor(ROOT.kRed)
    grx.Draw("p same")
    legx = TLegend(0.1, 0.7, 0.48, 0.9)
    legx.SetFillColorAlpha(0, 0.0)
    legx.SetLineColorAlpha(0, 0.0)

    c1.cd(4)
    gry0.SetMarkerStyle(20)
    gry0.SetMarkerSize(0.5)
    gry0.GetXaxis().SetTitle("Z/um")
    gry0.GetYaxis().SetTitle("Y/um")
    gry0.Draw("ap")
    gry.SetMarkerStyle(20)
    gry.SetMarkerSize(0.6)
    gry.SetMarkerColor(ROOT.kRed)
    gry.Draw("p same")
    legy = TLegend(0.1, 0.7, 0.48, 0.9)
    legy.SetFillColorAlpha(0, 0.0)
    legy.SetLineColorAlpha(0, 0.0)

    # Draw axes
    for i, params in enumerate(new_collections):
        X0, Y0, Z0, TX, TY = params
        fitXZ = TF1(f"fitXZ_{i}", "[0]+[1]*(x-[2])", zmin, zmax)
        fitXZ.SetParameters(X0, TX, Z0)
        fitXZ.SetLineColor(ROOT.kOrange + 9 if i == 0 else ROOT.kBlue)
        fitXZ.SetLineWidth(2 if i == 0 else 1)
        c1.cd(3)
        fitXZ.Draw("L SAME")
        if i == 0:
            legx.AddEntry(fitXZ, f"{labels[top_indices[i]]} (best)", "l")
        else:
            legx.AddEntry(fitXZ, labels[top_indices[i]], "l")

        fitYZ = TF1(f"fitYZ_{i}", "[0]+[1]*(x-[2])", zmin, zmax)
        fitYZ.SetParameters(Y0, TY, Z0)
        fitYZ.SetLineColor(ROOT.kOrange + 9 if i == 0 else ROOT.kBlue)
        fitYZ.SetLineWidth(2 if i == 0 else 1)
        c1.cd(4)
        fitYZ.Draw("L SAME")
        if i == 0:
            legy.AddEntry(fitYZ, f"{labels[top_indices[i]]} (best)", "l")
        else:
            legy.AddEntry(fitYZ, labels[top_indices[i]], "l")

    c1.cd(3)
    legx.DrawClone()

    c1.cd(4)
    legy.DrawClone()

    c1.SaveAs(t)