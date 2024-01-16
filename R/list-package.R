

#' The 1991 National Race and Politics Survey
#' 
#' This dataset is a subset of the 1991 National Race and Politics Survey and
#' contains the item count technique or the list experiment.  The main question
#' reads as follows: \verb{ Now I'm going to read you four things that
#' sometimes make people angry or upset.  After I read all (three/four), just
#' tell me HOW MANY of them upset you.  (I don't want to know which ones, just
#' how many.)  (1) "the federal government increasing the tax on gasoline;" (2)
#' "professional athletes getting million-dollar-plus salaries;" (3) "large
#' corporations polluting the environment;" (4) "black leaders asking the
#' government for affirmative action."}
#' 
#' where the last item is presented only with the treatment group and the
#' control list only contains the first three items.
#' 
#' 
#' @name affirm
#' @docType data
#' @format A data frame containing the following 6 variables for 1171
#' observations.
#' 
#' \tabular{llll}{ y \tab numeric \tab the number of items that make
#' respondents angry \tab 0 - 4 \cr south \tab numeric \tab whether or not a
#' respondents live in a southern state \tab 0 - 1 \cr male \tab numeric \tab
#' whether or not a respondent is male \tab 0 - 1 \cr college \tab numeric \tab
#' whether or not a respondent attended some college \tab 0 - 1 \cr age \tab
#' numeric \tab age of a respondent divided by 10 \cr treat \tab numeric \tab
#' treatment status \tab 0 - 1 }
#' @source The full data set is available at SDA (Survey Documentation and
#' Analysis; \url{https://sda.berkeley.edu/D3/Natlrace/Doc/nrac.htm})
#' @keywords datasets
NULL





#' Five List Experiments with Direct Questions
#' 
#' A dataset containing the five list experiments in Aronow, Coppock, Crawford,
#' and Green (2015)
#' 
#' 
#' @name combinedListExps
#' @docType data
#' @format A data frame with 1023 observations and 23 variables
#' @keywords datasets
NULL





#' The 2012 Mexico Elections Panel Study
#' 
#' This dataset is a subset of the 2012 Mexico Elections Panel Study and
#' contains a list experiment question. It reads as follows: \verb{ I am going
#' to read you a list of four activities that appear on this card and I want
#' you to tell me how many of these activities you have done in recent weeks.
#' Please don''t tell me which ones, just HOW MANY. The four activities are...
#' (SHOW CARD AND READ) a. See television news that mentions a candidate b.
#' Attend a campaign event c. Exchange your vote for a gift, favor, or access
#' to a service d. Talk about politics with other people}
#' 
#' where item c. is presented only to the treatment group, and the control list
#' only contains the other three items.
#' 
#' 
#' @name mexico
#' @docType data
#' @format A data frame containing the following variables for 1004
#' observations.  \tabular{llll}{ y \tab numeric \tab the number of items that
#' make respondents angry \tab 0 - 4 \cr
#' 
#' mex.t \tab numeric \tab treatment status \tab 0-1 \cr mex.male \tab numeric
#' \tab whether or not a respondent is male\tab 0-1 \cr mex.age \tab numeric
#' \tab age of a respondent \cr mex.education \tab numeric \tab respondent's
#' level of education \tab 0-9 \cr mex.y.all \tab numeric \tab the number of
#' activities that respondent did \tab 0 - 4 \cr mex.vote \tab numeric \tab
#' respondent's self-reported turnout \tab 0-1 \cr mex.age2 \tab numeric \tab
#' age of a respondent, squared \cr mex.interest \tab numeric \tab how
#' interested respondent is in politics \tab 1-4 \cr mex.married \tab numeric
#' \tab indicator for whether respondent is married \tab 0-1 \cr mex.pidpanw2
#' \tab numeric \tab indicator for whether respondent identifies with PAN party
#' \tab 0-1 \cr mex.pidprdw2 \tab numeric \tab indicator for whether respondent
#' identifies with PRD party \tab 0-1 \cr mex.pidpriw2 \tab numeric \tab
#' indicator for whether respondent identifies with PRI party \tab 0-1 \cr
#' mex.votecard \tab numeric \tab respondent's enumerator-verified turnout \tab
#' 0-1 \cr mex.urban \tab numeric \tab indicator for whether respondent lives
#' in urban area \tab 0-1 \cr mex.cleanelections \tab numeric \tab indicator
#' for whether respondent thinks elections were clean \tab 0-1 \cr
#' mex.cleanelectionsmiss \tab numeric \tab indicator for whether
#' cleanelections variable was missing \tab 0-1 \cr mex.metro \tab numeric \tab
#' indicator for whether respondent lives in Mexico City metro area \tab 0-1
#' \cr mex.centralregion \tab numeric \tab indicator for whether respondent
#' lives in Mexico's central region \tab 0-1 \cr mex.northregion \tab numeric
#' \tab indicator for whether respondent lives in Mexico's north region \tab
#' 0-1 \cr mex.wealth \tab numeric \tab scale for respondent's wealth, based on
#' household asset indicators \cr mex.epnapprove \tab numeric \tab respondent's
#' approval rating of Enrique Pena-Nieto \tab 1-11 \cr mex.havepropoganda \tab
#' numeric \tab indicator for whether respondent has propaganda outside their
#' home \tab 0-1 \cr mex.concurrent \tab numeric \tab indicator for whether
#' respondent lives in state with concurrent elections \tab 0-1 \cr mex.loyal
#' \tab numeric \tab indicator for whether respondent strongly identifies with
#' the PAN, PRI, or PRD party \tab 0-1 \cr mex.direct \tab numeric \tab
#' indicator for whether respondent directly reports an attempt to buy their
#' vote \tab 0-1 }
#' @source The full data set is available at the Mexico Panel Study website
#' (\url{http://mexicopanelstudy.mit.edu/})
#' @keywords datasets
NULL





#' The 1994 Multi-Investigator Survey
#' 
#' This dataset is a subset of the 1994 Multi-Investigator Survey and contains
#' the item count technique or the list experiment.  The main question reads as
#' follows: \verb{ Now I'm going to read you four things that sometimes make
#' people angry or upset.  After I read all (three/four), just tell me HOW MANY
#' of them upset you.  (I don't want to know which ones, just how many.)  (1)
#' "the federal government increasing the tax on gasoline;" (2) "professional
#' athletes getting million-dollar-plus salaries;" (3) "requiring seatbelts be
#' used when driving;" (4) "large corporations polluting the environment;" (5)
#' "black leaders asking the government for affirmative action."}
#' 
#' where the last item is presented only with the treatment group and the
#' control list only contains the first three items.
#' 
#' The survey also includes a question in which attitudes toward the sensitive
#' item are asked directly.
#' 
#' \verb{ Now I'm going to ask you about another thing that sometimes makes
#' people angry or upset. Do you get angry or upset when black leaders ask the
#' government for affirmative action?}
#' 
#' 
#' @name mis
#' @docType data
#' @format A data frame containing the following 6 variables for 1171
#' observations.
#' 
#' \tabular{llll}{ y \tab numeric \tab the number of items that make
#' respondents angry \tab 0 - 4 \cr sensitive \tab numeric \tab whether or not
#' the sensitive item (asked directly) makes respondents angry \tab 0 - 1 \cr
#' south \tab numeric \tab whether or not a respondents live in a southern
#' state \tab 0 - 1 \cr male \tab numeric \tab whether or not a respondent is
#' male \tab 0 - 1 \cr college \tab numeric \tab whether or not a respondent
#' attended some college \tab 0 - 1 \cr age \tab numeric \tab age of a
#' respondent divided by 10 \cr democrat \tab numeric \tab whether not a
#' respondent identifies as a Democrat \tab 0 - 1 \cr republican \tab numeric
#' \tab whether not a respondent identifies as a Republican \tab 0 - 1\cr
#' independent \tab numeric \tab whether not a respondent identifies as an
#' independent \tab 0 - 1 \cr treat \tab numeric \tab treatment status \tab 0 -
#' 1 \cr list.data \tab numeric \tab indicator for list experiment subset
#' (treatment and control groups) \tab 0 - 1 \cr sens.data \tab numeric \tab
#' indicator for direct sensitive item subset \tab 0 - 1 \cr
#' 
#' }
#' @source The full data set is available at SDA (Survey Documentation and
#' Analysis; \url{https://sda.berkeley.edu/D3/Multi/Doc/mult.htm})
#' @keywords datasets
NULL





#' The 1991 National Race and Politics Survey
#' 
#' This dataset is a subset of the 1991 National Race and Politics Survey and
#' contains a list experiment with two sensitive items.  The main questions
#' read as follows: \verb{ Now I'm going to read you four things that sometimes
#' make people angry or upset.  After I read all (three/four), just tell me HOW
#' MANY of them upset you.  (I don't want to know which ones, just how many.)
#' (1) "the federal government increasing the tax on gasoline;" (2)
#' "professional athletes getting million-dollar-plus salaries;" (3) "large
#' corporations polluting the environment;" (4) "a black family moving next
#' door to you."}
#' 
#' where the last item is presented only with the treatment group and the
#' control list only contains the first three items.
#' 
#' The second sensitive item replaces item (4) with \verb{ (4) "black leaders
#' asking the government for affirmative action." }
#' 
#' Treatment status one (treat == 1) is the "black family" item and status two
#' is the "affirmative action" item.
#' 
#' 
#' @name multi
#' @docType data
#' @format A data frame containing the following 6 variables for 1795
#' observations.
#' 
#' \tabular{llll}{ y \tab numeric \tab the number of items that make
#' respondents angry \tab 0 - 4 \cr south \tab numeric \tab whether or not a
#' respondents live in a southern state \tab 0 - 1 \cr male \tab numeric \tab
#' whether or not a respondent is male \tab 0 - 1 \cr college \tab numeric \tab
#' whether or not a respondent attended some college \tab 0 - 1 \cr age \tab
#' numeric \tab age of a respondent divided by 10 \cr treat \tab numeric \tab
#' treatment status \tab 0 - 2 }
#' @source The full data set is available at SDA (Survey Documentation and
#' Analysis; \url{https://sda.berkeley.edu/D3/Natlrace/Doc/nrac.htm})
#' @keywords datasets
NULL





#' The 1991 National Race and Politics Survey
#' 
#' This dataset is a subset of the 1991 National Race and Politics Survey and
#' contains the item count technique or the list experiment.  The main question
#' reads as follows: \verb{ Now I'm going to read you four things that
#' sometimes make people angry or upset.  After I read all (three/four), just
#' tell me HOW MANY of them upset you.  (I don't want to know which ones, just
#' how many.)  (1) "the federal government increasing the tax on gasoline;" (2)
#' "professional athletes getting million-dollar-plus salaries;" (3) "large
#' corporations polluting the environment;" (4) "a black family moving next
#' door to you."}
#' 
#' where the last item is presented only with the treatment group and the
#' control list only contains the first three items.
#' 
#' 
#' @name race
#' @docType data
#' @format A data frame containing the following 6 variables for 1213
#' observations.
#' 
#' \tabular{llll}{ y \tab numeric \tab the number of items that make
#' respondents angry \tab 0 - 4 \cr south \tab numeric \tab whether or not a
#' respondents live in a southern state \tab 0 - 1 \cr male \tab numeric \tab
#' whether or not a respondent is male \tab 0 - 1 \cr college \tab numeric \tab
#' whether or not a respondent attended some college \tab 0 - 1 \cr age \tab
#' numeric \tab age of a respondent divided by 10 \cr treat \tab numeric \tab
#' treatment status \tab 0 - 1 }
#' @source The full data set is available at SDA (Survey Documentation and
#' Analysis; \url{https://sda.berkeley.edu/D3/Natlrace/Doc/nrac.htm})
#' @keywords datasets
NULL



