.onAttach <- function(lib, pkg) { # .First.lib
	#library.dynam("galgo", pkg, lib)
    #dyn.load(paste("galgoDistance",.Platform$dynlib.ext,sep="")) 
	#lockEnvironment(as.environment("package:datasets"), TRUE)
	if(.Platform$OS.type == "windows" && interactive() && .Platform$GUI == "Rgui") addVigs2WinMenu("galgo")
	packageStartupMessage("galgo v1.2-01 (19-March-2014) was loaded.\n")
    packageStartupMessage("See 'packages' under R help for tutorial and manual.\n")
}

.onUnload <- function(libpath) { # .Last.lib = function(lib, pkg)
	#library.dynam.unload("galgo")
}

#THIS FUNCTION AS BEEN TAKEN AS IT IS FROM BIOBASE PACKAGE
addVigs2WinMenu <- function (pkgName) 
{
	vigs <- ""
    vigFile <- system.file(c("doc/Tutorial.pdf", "doc/Galgo.pdf"), package = pkgName)
    if (any(file.exists(vigFile))) {
		vigs <- vigFile[file.exists(vigFile)]
        #vigMtrx <- .readRDS(vigFile)
        #vigs <- file.path(.find.package(pkgName), "doc", vigMtrx[, "PDF"])
        #names(vigs) <- vigMtrx[, "Title"]
		names(vigs) <- c("Tutorial","Functions")[file.exists(vigFile)]
    }
    if (!"Vignettes" %in% winMenuNames()) 
        winMenuAdd("Vignettes")
    pkgMenu <- paste("Vignettes", pkgName, sep = "/")
    winMenuAdd(pkgMenu)
    for (v in 1:length(vigs)) {
		i <- vigs[v]
        item <- paste(names(vigs)[v],": ",basename(i),sep="") #sub(".pdf", "", basename(i))
        winMenuAddItem(pkgMenu, item, paste("shell.exec(\"", 
            as.character(i), "\")", sep = ""))
    }
}
