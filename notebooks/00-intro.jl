### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ 9802aaa7-3f9a-47b7-b6ab-90c4f39b7335
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ b62ba51a-b7a8-433c-84dd-bd7d221ffa3c
begin
	@quickactivate "ConvectionIsotopes"
	using DelimitedFiles
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("common.jl"))
	
	md"Loading modules for the ConvectionIsotopes project..."
end

# ╔═╡ fc7b6caa-6ced-11ec-0701-6f55729e22dc
md"
# 00. ColumbiaIsotope - Linking Isotopic data with Precipitation Measurements

In this project, we compare the isotopic station data for deuterium and 18-oxygen measurements with rainfall rates over various locations in Columbia.  However, things are moe easily said then done, and below is the outline:

1. Basic analysis of the overall average, and monthly climatology (both GPM and ERA5)
   * Filtered Land-Sea Mask for future use?  Allows us to split the GeoRegion into Atlantic, Pacific and Inland Regions?
   * Precipitation (both GPM and ERA5)
   * Vertical (and Horizontal??) Wind (ERA5 only)
2. Vertical profiles for different vertical-velocity-weighted mean pressure using ERA5 data
3. Basic overview of the station data including
   * Comparison between GPM and station rainfall (remember to adjust for timezones)
   * Looking at the Isotope data (there are errors in the file) and seeing if I can filter them out
   * Binned depletion ratio against precip, vertical-velocity-weighted pressure, is there enough data to conclude that 
4. Same as (3), but for WRF datasets
   * Binned depletion ratio against precip, vertical-velocity-weighted pressure
   * If (3) does not show much trends for depletion against vertical-velocity-weighted pressure, what is the length of observations needed?
"

# ╔═╡ 5ad5ac8f-f7fa-4a63-8600-ad93ae8096f6
begin
	coast = readdlm(datadir("GLB-i.txt"),comments=true,comment_char='#')
	x = coast[:,1]
	y = coast[:,2]
	md"Loading coastlines data ..."
end

# ╔═╡ 853a6546-4823-4cca-bf46-72cbc45c4df0
md"
### A. Compiling station information

We manually compiled a summary of the station location information in the project location `src/stninfo.csv`.  The format is .csv instead of .jld2 so that it can be read globally, as .nc doesn't really work with lists like this.

Stored information: (1) Name, (2) Longitude, (3) Latitude, (4) Height where available

There are also stations without longitude/latitude coordinates given, in which case I used a rough google-search to find the approximate location and put it in, so that we have a sense of where the stations are until we can (if possible) get the exact location.

Since San Andreas has both daily and monthly data, I've renamed the station in the daily data to be \"San Andres_01\" as with all the other daily data station names
"

# ╔═╡ 73ab7928-b428-4792-af0e-1f44b42beee9
begin
	infoall = stninfoall()
	infody  = stninfody()
	infomo  = stninfomo()
	md"Loading station location information ..."
end

# ╔═╡ 8bbbc366-114e-49a5-8762-8515a7f1de42
begin
	pplt.close(); f1,a1 = pplt.subplots(ncols=3,axwidth=2)
	
	a1[1].scatter(infoall[:,2],infoall[:,3])
	a1[2].scatter(infody[:,2],infody[:,3])
	a1[3].scatter(infomo[:,2],infomo[:,3])

	a1[1].format(ultitle="(a) All")
	a1[2].format(ultitle="(b) Daily")
	a1[3].format(ultitle="(c) Monthly")

	for ax in a1
		ax.plot(x,y,lw=0.5,c="k")
		ax.format(xlim=(-85,-70),ylim=(0,15),suptitle="Available Columbia Stations")
	end
	
	f1.savefig(plotsdir("00-stations.png"),transparent=false,dpi=300)
	load(plotsdir("00-stations.png"))
end

# ╔═╡ 03c0f4ae-e0d7-4ace-8cb5-cb015b8ebe2b
md"
### B. Splitting the Station Data

We have split the data given to us into two different xlsx files (we were given the files in this format), one with the summary of the station data, and the other with a more in-depth analysis of the monthly data (this one runs from 2013-2017 or thereabouts)

The station information is taken from the xlsx file (`data/IsotopeDataSummry.xlsx`) with the summary of the station data since it is most readily accessible and in the best tabular format and is thus the easiest to sort.  The file has two notebooks within, split into daily and monthly data.

It is to be noted that the San Andres station has two different sets of (I assume to be independently) collected data.  This is because this is the only station that has both monthly and daily datasets, but these datasets are collected over different time periods (daily data was collected later).  Also, the errors in the data collection in the daily data are also found in the San Andres data, which is why we believe that the daily and monthly datasets were collected independently even though the stations are basically in the same location.
"

# ╔═╡ Cell order:
# ╟─fc7b6caa-6ced-11ec-0701-6f55729e22dc
# ╟─9802aaa7-3f9a-47b7-b6ab-90c4f39b7335
# ╟─b62ba51a-b7a8-433c-84dd-bd7d221ffa3c
# ╟─5ad5ac8f-f7fa-4a63-8600-ad93ae8096f6
# ╟─853a6546-4823-4cca-bf46-72cbc45c4df0
# ╟─73ab7928-b428-4792-af0e-1f44b42beee9
# ╟─8bbbc366-114e-49a5-8762-8515a7f1de42
# ╟─03c0f4ae-e0d7-4ace-8cb5-cb015b8ebe2b
