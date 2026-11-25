"""Shared figure colors, typography, and publication exports."""
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.text import Text
from matplotlib.transforms import Affine2D, Bbox

PALETTE={
    "primary_dark":"#0B4F6C", "primary_teal":"#0D9488", "primary_light":"#5EEAD4",
    "risk_dark":"#B93815", "risk_orange":"#D95F02", "risk_soft":"#FDBA74",
    "anchor_purple":"#7E22CE", "anchor_lilac":"#C084FC", "neutral_dark":"#1E293B",
    "neutral_mid":"#64748B", "neutral_light":"#CBD5E1", "grid_bg":"#F1F5F9", "white":"#FFFFFF",
}
CMAP_DIV=LinearSegmentedColormap.from_list("jtm_div",["#0B4F6C","#38A3A5","#F8FAFC","#F4A261","#B93815"])
CMAP_SEQ=LinearSegmentedColormap.from_list("jtm_seq",["#F1F5F9","#99F6E4","#0D9488","#0B4F6C","#082F49"])
CMAP_RISK=LinearSegmentedColormap.from_list("jtm_risk",["#F8FAFC","#FDE68A","#F97316","#B93815","#7F1D1D"])
CMAP_UNITS=LinearSegmentedColormap.from_list("sampling_units",[
    "#0B4F6C","#0D9488","#5EEAD4","#F0F9F6","#FDBA74","#D95F02","#B93815",
])


def apply_jtm_theme():
    mpl.rcParams.update({
        "figure.facecolor":"white", "axes.facecolor":"white", "axes.edgecolor":"#334155",
        "axes.linewidth":.75, "axes.spines.top":False, "axes.spines.right":False,
        "font.family":"sans-serif", "font.sans-serif":["Arial","Helvetica","DejaVu Sans"],
        "text.color":PALETTE["neutral_dark"], "axes.labelcolor":PALETTE["neutral_dark"],
        "xtick.color":PALETTE["neutral_dark"], "ytick.color":PALETTE["neutral_dark"],
        "font.size":7.5, "axes.labelsize":7.5, "xtick.labelsize":6.8,"ytick.labelsize":6.8,
        "legend.fontsize":6.8,"legend.frameon":False,"pdf.fonttype":42,"ps.fonttype":42,
        "svg.fonttype":"none","savefig.dpi":300,"axes.unicode_minus":True,
    })


def inset_content(fig,fraction=.01):
    """Move the complete layout inward without changing the page or type sizes."""
    fx,fy=fraction if isinstance(fraction,tuple) else (fraction,fraction)
    placement=Affine2D().scale(1-2*fx,1-2*fy).translate(fx,fy)
    for axes in fig.axes:
        bounds=axes.get_position(original=True)
        axes.set_position(Bbox(placement.transform(bounds.get_points())))
    for artist in fig.texts:
        artist.set_position(placement.transform(artist.get_position()))
    for artist in fig.artists:
        artist.set_transform(placement+fig.transFigure)


def export(fig,stem):
    # Instantiate tick labels before enforcing the minimum type size.
    fig.canvas.draw()
    for item in fig.findobj(Text):
        if item.get_text():
            item.set_fontsize(max(6.5,item.get_fontsize()))
    for extension in (".pdf",".svg",".png"):
        fig.savefig(stem.with_suffix(extension),dpi=300,facecolor="white",bbox_inches=None,pad_inches=0)
