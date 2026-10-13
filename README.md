# dssd Survey Designer

A browser-based (WebAssembly / Shinylive) front end for the
[dssd](https://github.com/DistanceDevelopment/dssd) package for distance
sampling survey design (Marshall & Rexstad).

Workflow: region (upload shapefile/GeoPackage/GeoJSON or use a dssd example)
→ design (line or point transects) → generate transects → coverage simulation
→ export (shapefile, GPX, CSV, coverage statistics, reproducible R script).

Everything runs in the user's browser; uploaded files are not sent to a server.
The first load takes a while because R and the packages are downloaded into
the browser.

## Run locally in R

```r
install.packages(c("shiny", "bslib", "DT", "sf", "dssd", "ggplot2", "zip"))
shiny::runApp()
```

## Build the WebAssembly version locally

```r
install.packages(c("shinylive", "httpuv"))
shinylive::export(".", "site")
httpuv::runStaticServer("site")
```

## Deployment

Pushing to `main` runs `.github/workflows/deploy.yml`, which exports the app
with shinylive and publishes it to GitHub Pages.
