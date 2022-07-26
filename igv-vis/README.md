## IGV Distribution

## Summary
Modified version of IGV for visualzing m6dA footprints on single molecules / single reads. m6dA footprints must be encoded in ML and MM BAM tags using the specification "N+a", and can be visualized by right-clicking on an alignment track and selecting the option "Color by" --> "Base Modification"

## Installation
IGV requires Java 11 in order to run. Java 11 can be installed using any package manager. 

### Installing Java 11
#### macOS:
Install <code>openjdk11</code> using [homebrew](https://formulae.brew.sh/formula/openjdk@11):

    brew install openjdk@11

Then link the newly installed jdk installation to the default JDK location on macOS

on M1 macOS:

    ln -s /opt/homebrew/Cellar/openjdk@11/11.0.15/libexec/openjdk.jdk /Library/Java/JavaVirtualMachines/openjdk.jdk

on x86 macOS: 

    ln -s /usr/local/Cellar/openjdk@11/11.0.15/libexec/openjdk.jdk /Library/Java/JavaVirtualMachines/openjdk.jdk

#### Linux
Install  <code>openjdk11</code> using <code>apt-get</code>:

    sudo apt-get install openjdk-11-jdk


## Running IGV
To launch IGV, source or run the script appropriate for your operating system:

* macOS - [igv.command](./IGV-dist/igv.command)
* Linux - [igv.sh](./IGV-dist/igv.sh)
* Windows - [igv.bat](./IGV-dist/igv.bat)


## Example Visualizations
A visualization of m6dA footprints on SMRT-Tag reads derived from the MYC-amplified osteosarcoma cell line [OS152](https://doi.org/10.1158/2159-8290.CD-17-1152) is provided in the folder [reasd](./reads/).

To load the visualization, open the session file [OS152_MYC.xml](./OS152_MYC.xml) with the modified IGV.


## Software versions
* [IGV](https://github.com/igvteam/igv/) at commit [5705dce6b36f88cf1a252ee158d9d96648c5a20e](https://github.com/igvteam/igv/commit/5705dce6b36f88cf1a252ee158d9d96648c5a20e)