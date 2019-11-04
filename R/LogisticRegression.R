LogisticRegression <- setRefClass("LogisticRegression", contains="Classifier",
  
  methods = list(
		initialize = function(...) {
			callSuper(...)
      buildClassifier()
      .self
		},

		buildClassifier = function() {
			.self$model <- glm(y~., family=binomial,
        data=data.frame("y"=factor(.self$y), t(.self$x), check.names=FALSE),
        maxit=50)
		},
    
    subselect = function(method="stepAIC") {
      switch(method,
             "stepAIC"=reduce.stepaic(),
             "anneal"=reduce.anneal(),
             "lasso"=reduce.lasso())
    },
    
    reduce.stepaic = function() {
      .self$model <- stepAIC(.self$model, trace=0, direction="backward")
    },
    
    reduce.anneal = function() {
      matr <- glmHmat(.self$model)
      .self$model <- anneal(matr$mat, kmin=3, H=matr$H, r=matr$r)
    },
    
    reduce.lasso = function() {
      .self$model <- glmnet(t(.self$x), .self$y, family="binomial", alpha=1)
    },
  
  reduce.mclust.stepaic = function() {
    mcl$new()
  },

		trainingPrediction = function() {
			.self$predicted.value <- predict.glm(.self$model, newdata=data.frame(t(.self$x),check.names=F), type="response")
			levels(.self$predicted.class) <- levels(.self$y)

			if(!is.null(.self$optimalCutoff)) {
				.self$predicted.class[.self$predicted.value > .self$optimalCutoff] <- levels(.self$y)[2]
			} else {
				.self$predicted.class[.self$predicted.value > 0.5] <- levels(.self$y)[2]
			}
		},

		testPrediction = function() {
			.self$predicted.value <- predict.glm(.self$model, newdata=data.frame(t(.self$x.test),check.names=F), type="response")
			levels(.self$predicted.class) <- levels(.self$y.test)

			if(!is.null(.self$optimalCutoff)) {
				.self$predicted.class[.self$predicted.value > .self$optimalCutoff] <- levels(.self$y.test)[2]
			} else {
				.self$predicted.class[.self$predicted.value > 0.5] <- levels(.self$y)[2]
			}
		},


		trainingStatistics = function() {
			.self$trainingPrediction()
			.self$computeROC()
			.self$computeAUC()
			.self$computeoptimalCutoff()
			.self$computeSensitivity()
			.self$computeSpecificity()
		},
    
		testStatistics = function() {
		  .self$testPrediction()
		  .self$computeROC()
		  .self$computeAUC()
		  .self$computeoptimalCutoff()
		  .self$computeSensitivity()
		  .self$computeSpecificity()
		}
		))
