# ygl

## Description

This script extracts and plots subcellular localization analysis data from [\[1\]](https://www.pnas.org/content/110/39/15842).
The data are retreived from this [server](http://128.179.34.6/twiki/bin/view/CellImaging/WebHome). There is also a [README](http://128.179.34.6/MMS_screen/datafiles/README.txt) explaining the various parameters found in the `.m` data. The plots are in interactive form and can be found [here](https://acubesat.gitlab.io/su/biology/yeast-biology-pages/dashboard.html).

## Data

For each strain, a single-cell rectangular image was captured across a total of 40 timepoints. To distinguish between very fine localization patterns as usually defined in cell biology and more objective geometrical shapes, the authors analyzed the images using a custom-made propbabilistic classification scheme. Six different shapes/patterns were defined:

1. **Periphery**: A fine outline of the cell contour, generally very well distinguishable. Representative strains include membrane proteins uniformly associated with the cell membrane, or in some cases bright dots distributed on the membrane.
2. **Structure**: This shape includes filaments, circles, and shape-forming dots that are often a direct indication for organelle-related localization of the protein.
3. **Punctate**: A number of distinct small dots of sizes smaller than 1μm (< 20% of the size of the cells). Typical representatives are *actin*, *lipid particles* and *peroxisomes*.
4. **Disk**: One dominant area of GFP signal contained in the interior of the cell. The diameter of these objects is at least around 25% of the diameter of the cell. Typical representatives of this group are strains localized in the *nucleus* and *nucleolus*, but also proteins in the *vacuole* or *vacuolar membrane*.
5. **Corona**: Broad ring (donut) around the center that can also be more sickle-shaped. Typical localizations that have a corona-like appearance are *cytoplasm* and in some cases *ER*.
6. **Homogeneous**: Cells where the fluoresence is uniformly distributed. In many cases homogeneous cells are of low intensity reflecting background levels

Because the boundaries between these shapes can be funny, the probabilistic classifier scheme is used to assign to each cell a probability vector reflecting the likelihood to belong to each corresponding pattern. For more info, visit the [SI](http://128.179.34.6/twiki/pub/CellImaging/SuppMaterial/Denervaud_Supplement.pdf).

---

[1] Dénervaud, N., Becker, J., Delgado-Gonzalo, R., Damay, P., Rajkumar, A. S., Unser, M., ... & Maerkl, S. J. (2013). A chemostat array enables the spatio-temporal analysis of the yeast proteome. Proceedings of the National Academy of Sciences, 110(39), 15842-15847.
