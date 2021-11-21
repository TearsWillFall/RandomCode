usethis::create_package(paste0(getwd(),"/RandomCode"))

)
setwd(paste0(getwd(),"/RandomCode/tools"))
setwd(paste0(getwd(),"/RandomCode"))
devtools::load_all()
devtools::document()


system("wget https://github.com/Illumina/strelka/archive/refs/tags/v2.9.10.tar.gz")

.libPaths()[1]
"manta/bin/configManta.py"
