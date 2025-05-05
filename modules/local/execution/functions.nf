def colored_outputs() {
    ANSI_RESET = "\u001B[0m";
    ANSI_BLACK = "\u001B[30m";
    ANSI_RED = "\u001B[31m";
    ANSI_GREEN = "\u001B[32m";
    ANSI_YELLOW = "\u001B[33m";
    ANSI_BLUE = "\u001B[34m";
    ANSI_PURPLE = "\u001B[35m";
    ANSI_CYAN = "\u001B[36m";
    ANSI_WHITE = "\u001B[37m";
    ANSI_BOLD = "\u001B[1m";

    def print_red = {  str -> ANSI_RED + str + ANSI_RESET }
    def print_black = {  str -> ANSI_BLACK + str + ANSI_RESET }
    def print_green = {  str -> ANSI_GREEN + str + ANSI_RESET }
    def print_yellow = {  str -> ANSI_YELLOW + str + ANSI_RESET }
    def print_blue = {  str -> ANSI_BLUE + str + ANSI_RESET }
    def print_cyan = {  str -> ANSI_CYAN + str + ANSI_RESET }
    def print_purple = {  str -> ANSI_PURPLE + str + ANSI_RESET }
    def print_white = {  str -> ANSI_WHITE + str + ANSI_RESET }
    def print_bold = { str -> ANSI_BOLD + str + ANSI_RESET }
}

def processChannelStatistics(ch_statistic) {
    return ch_statistic
        .groupTuple(by: 0)
        .map { sample, channels, counts ->
            def mappedCounts = [:] // Create an empty map to hold channel:count mappings
            channels.eachWithIndex { channel, i ->
                mappedCounts[channel] = counts[i] // Map each channel to its corresponding count
            }
            return [sample, mappedCounts]
        }
        .toList()
        .map { allSamples ->
            def maxWidths = allSamples.collect { it[0].toString().length() }.max()
            def channelWidths = allSamples*.get(1).collectMany { it.keySet() }.unique().collectEntries { [(it): it.toString().length()] }
            allSamples.each { sample, counts -> counts.each { k, v -> channelWidths[k] = Math.max(channelWidths[k], v.toString().length()) } }
            def header = "sample".padRight(maxWidths) + "\t" + channelWidths.collect { k, v -> k.padRight(v) }.join("\t")
            def rows = allSamples.collect { sample, counts ->
                def row = sample.toString().padRight(maxWidths) + "\t" + channelWidths.collect { k, v -> counts.get(k, "0").toString().padRight(v) }.join("\t")
                return row
            }
            return ([header] + rows).join("\n")
        }
        // .set { sample_statistic_table }
}

def processMergedStatisticsChannel(ch_statistic_merged) {
    return ch_statistic_merged
        .groupTuple(by: 0)
        .map { sample, channels, counts ->
            def mappedCounts = [:]                      
            channels.eachWithIndex { channel, i ->
                mappedCounts[channel] = counts[i]       
            }
            def stats = mappedCounts.collect { k, v -> "$k: $v" }.join(", ")
            return "$sample: $stats"
        }
}