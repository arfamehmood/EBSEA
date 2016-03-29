filterGenes <- function(x, fc=1.25, fdr=0.01)
{

    # Checking Parameters
    message("Performing Filtering")
    if(is.null(x)){stop('There is no gene data....')}

    FC <- apply(matrix(colnames(x)), 1, function(x)if(x == 'FC') FC='found')

    if(is.null(FC)){stop('There are no Fold Changes for data')}
    FDR <- apply(matrix(colnames(x)), 1, function(x)if(x == 'FDR') FDR='found')
    if(is.null(FDR)){stop('There are no FDR values in the data')}

    # Obtaining Up and down regulated genes
    urGenes <- x[x$FC > fc & x$FDR < fdr, ]
    drGenes <- x[x$FC < -fc & x$FDR < fdr, ]
    return(list("Upregulated" = urGenes[with(urGenes, order(-FC)), ],
                "Downregulated" = drGenes[with(drGenes, order(-FC)), ]))

}
