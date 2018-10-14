\name{ALL}
\alias{ALL}
\docType{data}
\title{ Acute Lymphoblastic Leukemia data (Yeoh et. al., 2002) for GALGO package }
\description{
  Acute Lymphoblastic Leukemia for GALGO package data published by Yeoh et. al.
  The original 360 pediatric acute leukemia samples were filtered by class. 233 samples are included corresponding to 5 classes EMLLA, HYP+50, MLL, T, and TEL.
  The Affymetrix microarray HG_U95Av2 containing 12,600 probesets were filteres by range and standard deviation resulting in 2,435 probesets (genes).
}
\usage{data(ALL)}
\details{ ALL data is complmented by ALL.classes which contain the classes for each column sample. }
\format{
  The format is:
         2,435 Rows : Genes
        233 Columns : Samples
          Row Names : ProbeId
	  Col Names : Samples Id
	  5 Classes : Lukemia Types: "EMLLA"=E2A-PBX1, "T"=T-ALL, "HYP+50"=Hyperdiploid > 50 Chromosomes, "MLL"=MLL rearragment, and "TEL"=TEL-AML1
}
\references{
  Eng-Juh Yeoh, Mary E. Ross, Sheila A. Shurtleff, W. Kent Williams, Divyen Patel, Rami Mahfouz, Fred G. Behm, Susana C. Raimondi, Mary V. Relling, Anami Patel, Cheng Cheng, Dario Campana,, Dawn Wilkins, Xiaodong Zhou, Jinyan Li, Huiqing Liu, Ching-Hon Pui, William E. Evans, Clayton Naeve, Limsoon Wong, and James R. Downing. \emph{Classification, subtype discovery, and prediction of outcome in pediatric acute lymphoblastic leukemia by gene expression profiling}. Cancer Cell. March 2002.
}
\examples{
data(ALL)
data(ALL.classes)
}
\keyword{datasets}


