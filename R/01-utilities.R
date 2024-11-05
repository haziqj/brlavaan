get_lav_stuff <- function(fit) {
  # Utility function to extract lavaan stuff
  list(
    lavmodel       = fit@Model,
    lavsamplestats = fit@SampleStats,
    lavdata        = fit@Data,
    lavoptions     = fit@Options
  )
}
