%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  
% 
% on Wed Feb 08 14:38:13 2006.
% 
% Generator was the Rdoc class, which is part of the R.oo package written
% by Henrik Bengtsson, 2001-2004.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{predict.BigBang}
\alias{predict.BigBang}
\alias{BigBang.predict}
\alias{predict.BigBang}
\alias{predict,BigBang-method}
\alias{predict}
\keyword{methods}
\keyword{internal}


 \title{Predicts the class or fitting of new set of samples}

 \description{
  Predicts the class or fitting of new set of samples.
 }

 \usage{\method{predict}{BigBang}(object,
	newdata,
	permanent=FALSE,
	newClass="UNKNOWN",
	func=classPredictionMatrix,
	scale=if (is.null(o$data$scale)) FALSE else o$data$scale,
	...)}

 \arguments{
	\item{newdata}{Matrix or data frame with the same dimensions than original data. columns are new samples, rows are genes.}
	\item{permanent}{Should the newdata become permenent in the bigbang object. This is used to make plots including the new samples. However, the new samples are added supposing a new class "UNKNOWN" which may be annoying because it is a false class (added to simplify the implementation).}
	\item{newClass}{Name of the new class or names of each new sample if they want to be distinguished in plots. It could be a numeric vector, which is interpreted as levels of original classes. The default is "UNKNOWN".}
	\item{func}{The function wanted to be called. If permament is TRUE, this is unnecessary. The default is classPredictionMatrix, which will produce the class prediction for the new data. If the BigBang object is not a classification problem this function is the desired fitting/predicting function. If it is not a function, no call be performed.}
	\item{...}{Further parameters passed to func.}
 }

 \value{
	Returns the result of the call to func.
 }

 \examples{\dontrun{
   #bb is a BigBang object
   #nd is a the new data frame, rows=genes, cols=samples
   cpm <- predict(bb, newdata=nd)
   cpm
 
   #permanent data = PLOTS
   cpm <- predict(bb, newdata=nd, permanent=TRUE)
   plot(bb, cpm, type="confusion")
   }
 }

 \references{Goldberg, David E. 1989 \emph{Genetic Algorithms in Search, Optimization and Machine Learning}. Addison-Wesley Pub. Co. ISBN: 0201157675}

 \author{Victor Trevino. Francesco Falciani Group. University of Birmingham, U.K. http://www.bip.bham.ac.uk/bioinf}

 \seealso{
   For more information see \code{\link{BigBang}}.
 }

 \keyword{methods}
