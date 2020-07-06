# misaem 1.0.0

A major release mainly update the follows:

* Add linear regression with missing values, see the function `miss.lm`.

* The main function for logistic regression is changed to `miss.glm`.

* The functions mimics the structure of widely used functions `lm` and `glm` for the case without missing values.

* Update misaem tutorial.

* Update README.md.

* Typos fixed.

# misaem 0.9.2

A minor release mainly fixing bug:

* Replace `if (class(obj) == "data.frame")  { ..... }` with `if (is(obj, "data.frame"))  { ..... }`. As within the formal S4 class system, a matrix `obj` could well be of a class that extends matrix, such as "matrix" "array" (length = 2).

# misaem 0.9.1

A minor release mainly fixing bugs and typos:

* Fix a bug in `model_selection`, now it can correctly perform model selection if the full model is the best model.

* In `pred_saem`, two methods for prediction of test set with missingness are provided.

* Fix some typos in the vignettes. The length of codes now fits the page wide of html.

* Delete unused Imports `magrittr` in DESCRIPTION file.

* Change the index of vignettes from 'SAEM' to 'misaem tutorial'.

* Update README.md.
