##
## This file contains the matplotlib style settings for generating
## noise requirement figures with normal colors.
##
## For a description of the settings and list of available options see:
## https://matplotlib.org/users/customizing.html#a-sample-matplotlibrc-file
##

## Figure canvas settings ##
figure.facecolor:       white
figure.edgecolor:       white
figure.figsize:         12, 9

## Plot axes settings ##
axes.labelsize:         8
axes.facecolor:         white
axes.edgecolor:         black
axes.labelcolor:        black
axes.labelweight:       normal
axes.labelpad:          6
axes.linewidth:         1
axes.grid:              True
axes.grid.axis:         both
axes.grid.which:        both
axes.titlesize:         8

## Text settings ##
font.size:              8
text.color:             black
text.usetex:            False
# [The below settings have effect only if usetex=False]
font.weight:            normal
font.family:            serif
font.serif:             Georgia
mathtext.fontset:       dejavuserif

## X-tick settings ##
xtick.labelsize:        8
xtick.color:            black
xtick.top:              True
xtick.bottom:           True
xtick.direction:        in
xtick.major.size:       4
xtick.major.width:      0.5
xtick.major.pad:        5
xtick.minor.size:       2
xtick.minor.width:      0.5
xtick.minor.visible:    True

## Y-tick settings ##
ytick.labelsize:        8
ytick.color:            black
ytick.left:             True
ytick.right:            True
ytick.direction:        in
ytick.major.size:       4
ytick.major.width:      0.5
ytick.major.pad:        3.5
ytick.minor.size:       2
ytick.minor.width:      0.5
ytick.minor.visible:    True

## Grid line settings ##
grid.color:             xkcd:cement
grid.linestyle:         -
grid.linewidth:         0.5
grid.alpha:             0.3

## Plot line settings ##
lines.linewidth:        2

## Legend settings ##
legend.borderpad:       0.2
legend.fancybox:        True
legend.fontsize:        8
legend.framealpha:      0.8
legend.handletextpad:   0.5
legend.labelspacing:    0.33
legend.loc:             best

patch.force_edgecolor:  True

## Save settings ##
savefig.dpi:            140
savefig.bbox:           tight
pdf.compression:        9
