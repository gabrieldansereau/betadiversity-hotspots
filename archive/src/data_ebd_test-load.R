#### Test exported EBD data ####
# !!! CAREFUL !!! Some steps might overflow memory

source(file.path("src", "required.R"))

# Load data
warblers <- read_tsv(here("..", "data", "ebd", "processed", "ebd_warblers_cut.csv"))
names(warblers) = gsub(" ", ".", names(warblers))
head(warblers)
warblers_sampling <- read_tsv(here("..", "data", "ebd", "processed", "ebd_warblers_sampling.csv"))
names(warblers_sampling) = gsub(" ", ".", names(warblers_sampling))
head(warblers_sampling)

# Check summary
# warblers_summary <- summary(warblers) # takes some time
# warblers_summary

# Check categories
warblers %>%
  filter(CATEGORY == "form") %>%
  select(CATEGORY, SCIENTIFIC.NAME, SUBSPECIES.SCIENTIFIC.NAME) %>%
  head(50)
warblers %>%
  filter(CATEGORY == "intergrade") %>%
  select(CATEGORY, SCIENTIFIC.NAME, SUBSPECIES.SCIENTIFIC.NAME) %>%
  head(50)
warblers %>%
  filter(CATEGORY == "issf") %>%
  select(CATEGORY, SCIENTIFIC.NAME, SUBSPECIES.SCIENTIFIC.NAME) %>%
  head(50)

# Check species
unique(warblers$SCIENTIFIC.NAME) # 63
# Check subspecies
unique(warblers$SUBSPECIES.SCIENTIFIC.NAME) # 41
# Check counts
sort(unique(warblers$OBSERVATION.COUNT)) # X is weird
warblers %>%
  filter(OBSERVATION.COUNT != "X") %>%
  pull(OBSERVATION.COUNT) %>%
  as.numeric %>%
  hist

# Check longitude
filter(warblers, LONGITUDE > 0) # possibly an error

warblers_strict <- warblers %>%
  filter(CATEGORY == "species",
       OBSERVATION.COUNT != "X",
       LONGITUDE < 0,
       PROTOCOL.TYPE == "Traveling",
       APPROVED == 1) %>%
  droplevels
summary(warblers_strict)

# Transform dates
tmp <- head(warblers, 10)
tmp %>%
  mutate_at(vars(OBSERVATION.DATE), list(YEAR = year, MONTH = month, DAY = day))
warblers <- warblers %>%
  mutate_at(vars(OBSERVATION.DATE), list(YEAR = year, MONTH = month, DAY = day))

# Check dates
warblers %>%
  count(YEAR) %>%
  print.data.frame()
warblers %>%
  count(MONTH) %>%
  print.data.frame()
warblers %>%
  count(DAY) %>%
  print.data.frame()
