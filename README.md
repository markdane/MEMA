# MEMA

The MEMA package includes functions to preprocess, QA, normalize and explore Microenvironment Microarray (MEMA) experiments. MEMAs are used in Oregon Health and Science's (OHSU) MEP-LINCS project. An over veiw of the project along with all data and protocols are avilable at https://www.synapse.org/#!Synapse:syn2862345/wiki/72486.

<br>

####Relative File Structure

The functions in this package assume the data is stored on the MEP-LINCS Synapse website or in a relative file structure that mimics the original server structure as described below:

* top level directories named by each well plate's barcode value  

  * A subdirectory named "Analysis"  
     an2omero_barcode.tsv   contains metadata about the cell lines and reagents in the well plate
     barcode_level1.tsv  contains cell level data and metadata  
     barcode_Lebel2.tsv  contains spot level data and metadata
       
  *  A subdirectory under "Analysis" named by the segmentation pipeline version such as "v2"  
      Raw data files from the segmentation pipeline  
      
* a directory named "study"  
  * a subdirectory named by a study name such as "HMEC122L_SS1"  
    studyName_Level3.tsv  spot level data from multiple plates that has been normalized  
    studyName_Level4.tsv replicate level data for a study that has been summarized from the level 3 data  
    
<br>
  
####Preprocessing Pipeline Overview
The functions in this package support data from both CellProfiler and GE InCell segmentation pipelines. The metadata can come OHSU's An! metadata database or excel files. MEMAs can be in 8 or 96 well plates. 

<br>

####Package Downloads
-   The MEMA package can be downloaded and installed from this repo with the command:

    ``` r
    devtools::install_github("MEP-LINCS/MEMA")

    ```
or the development version with the command:

    ``` r
    devtools::install_github("MEP-LINCS/MEMA", ref="develop")

    ```