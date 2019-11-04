RandomForest <- setRefClass("RandomForest", contains="Classifier",
  
  methods = list(
		initialize = function(...) {
			callSuper(...)
      buildClassifier()
      .self
		},

		buildClassifier = function() {
			yx <- data.frame("y"=factor(.self$y), t(.self$x))
			colnames(yx)[2:ncol(yx)] <- rownames(.self$x)
			.self$model <- randomForest(y~., data=yx)
		},
    
		trainingPrediction = function() {
			.self$predicted.value <- predict(.self$model, newdata=data.frame(t(.self$x),check.names=F), type="prob")
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
