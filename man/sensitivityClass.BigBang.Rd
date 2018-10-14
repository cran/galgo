%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:08 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{sensitivityClass.BigBang}
\alias{sensitivityClass.BigBang}
\alias{BigBang.sensitivityClass}
\alias{sensitivityClass.BigBang}
\alias{sensitivityClass,BigBang-method}
\alias{sensitivityClass}
\keyword{methods}
\keyword{internal}


 \title{Computes the sensitivity of class prediction}

 \description{
  Computes the sensitivity of class prediction.
 }

 \usage{\method{sensitivityClass}{BigBang}(o, cm, ...)}

 \arguments{
	\item{cm}{The confusion matrix or the class prediction matrix. If missing, \code{confusionMatrix} method is called using the object and \code{...} as other arguments}
	\item{..}{Further parameters when \code{cm} is missing.}
 }

 \details{
  Sensitivity is the probability that a sample of class \code{X} will be predicted as the same class \code{X}. High sensitivity detect true positives.
  Sensitivity = TP / (TP + FN)
  TP - True Positives: Example for class A, TP = Paa
  FN - False Negatives: Example for class A, FN = Pab + Pac + Pax
  Confusion Matrix:
				[ Predicted Class ]
			ClassA	ClassB	ClassC	"misclass"
  ClassA	Paa		Pab		Pac		Pax
  ClassB	Pba     Pbb		Pbc		Pbx
  ClassC	Pca     Pcb		Pcc		Pcx
 }

 \value{
  A vector with the sensitivities of prediction for every class.
 }

 \examples{\dontrun{
   #bb is a BigBang object
   cpm <- classPredictionMatrix(bb)
   cpm
   cm <- confusionMatrix(bb)
   cm
   #equivalent and quicker because classPredictionMatrix is provided
   cm <- confusionMatrix(bb, cpm)
   cm
 
   specificityClass(bb, cm)
   specificityClass(bb, cpm)
   specificityClass(bb)
   # all are equivalent
   sensitivityClass(bb, cpm)
   sensitivityClass(bb, cm)
   sensitivityClass(bb)
   # all are equivalent
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
   \code{\link[galgo:classPredictionMatrix.BigBang]{*classPredictionMatrix}()},
   \code{\link[galgo:confusionMatrix.BigBang]{*confusionMatrix}()}.
 }

 \keyword{methods}
