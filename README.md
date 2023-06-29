# Cluster-Fitting
MATLAB code for fitting clusters in segmented cells. 
Readme file for fitting clusters.
General: This script works by manually drawing cell shapes, and the fitting clusters within each cell. See methods section of “Visual detection of submicroscopic protein clusters with a phase-separation-based fluorescent reporter” for more information. 
Segment_and_fit_clusters.m: This is the main script that can be used to fit cells. It has several sections. 
•	Drawing cells
o	Autoload files enables importing of files by constructing the filename(s) in an image set and looping over them to import all desired files
o	Load files (manual) opens a file browser window where files can manually be selected to import them.  
o	Pick image, set lookup table displays an image (image number can be changed to select which image to display) and the lookup table can be set using variables imin and imax. 
o	Draw a new cell enables the user to draw a cell, selecting a grouping of pixels for fitting clusters within. The cell will not be saved until the keep chunk of code is run after a cell is drawn. 
o	Redraw will zoom the image in on the current cell that is being drawn, enabling more precise drawing of the cell. 
o	Keep saves a drawn cell. Cell will be outlined in red after it is drawn, but outlined in green once it is saved. 
o	Modify a cell enables selection and redrawing of a cell that has been saved while preserving it’s cell number.
o	Delete a cell enables selection of a previously saved cell to be deleted.
o	Pass fits from image x to image y will copy cells saved on one image to another image. 
o	Import fits from file will import cell fits previously saved into the current workspace. 
o	Save fits to file will save the current cell fits (on all images) to a file where they can later be imported. 
o	Notes on variable names: fits for the entire image set are saved into variable ‘cell_inds_master’ which is a cell array with one entry for each image in the loaded stack. Each entry will contain another cell array with an entry for each cell. Each entry in these nested cell arrays will be a 1xn matrix containing the linear indices of each pixel in the drawn cell (where n is the number of pixels in the drawn cell). The variable ‘cell_inds’ is this nested cell array (one entry for each drawn cell) for the currently opened image. 
•	Fitting clusters
o	Clusters and centers on current image: This chunk of code will fit clusters within cells drawn on the currently displayed image. The main purpose of this is to find appropriate values for the variables “center_parameter” and “cluster_parameter”. Center parameter determines how bright a cluster has to be to be fit, and cluster parameter determines how bright a pixel must be to be included in a cluster. 
o	Batch clusters and centers: This code will fit clusters in every cell in every image in the stack. 
o	Cluster and outline display: This code can be used to display images with fit clusters after batch fitting. 
o	Data export and masks: This code will export a “cell” spreadsheet, containing x position, y position, intensity, cell of origin, and image of origin for each pixel included in a cell. A “cluster” spreadsheet will also be exported, containing the x positon, y position, intensity, cell of origin, image of origin, and cluster of origin for each pixel fit as a cluster across the image set. An “info” spreadsheet will also be exported, linking “image_number” to the filename of this image, as well as recording the parameters used to fit the clusters. 
