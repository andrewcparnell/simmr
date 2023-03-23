# Check global variables
if (getRversion() >= "2.15.1") {
  utils::globalVariables(
    names = c(
      "Proportion",
      "Source",
      "Type",
      "..level..",
      "density"
    ),
    package = "simmr",
    add = FALSE
  )
}
