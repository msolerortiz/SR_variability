
######### DATA CONFIG ####

# configuration variables
G_sample_frequency = 187
G_data_file_path = "/media/msolerortiz/ExternalDrive/"
G_debug = FALSE

######### FUNCTIONS ####

# ----- data structure -----
# Data files are Matlab files stored in an external drive with a tree structure 
# YYYY/MM and the following format: 
#
# c_YYYY_MM_DD_HHMMSS_OR_DES.mat
#
# with YYYY_MM_DD_HHMMSS being the timestamp of the data, ORientation being 
# either EW or NS and DEScriptor indicating preprocessing status.
#

get_datetime_from_file <- function(file) {
  raw_info = unlist( strsplit(file, "_") )
  file_date = paste(raw_info[2],raw_info[3],raw_info[4], sep="/")
  file_datetime = paste(file_date, raw_info[5], sep=" ")
  #file_name is arrival date, 30 minutes shifting to get start date.
  date_info = ymd_hms(file_datetime) - minutes(30)
  return(date_info)
}

char_and_add_zero <-function(value) {
  if (value > 9) {
    value_as_char = as.character(value)
  } 
  
  else {
    value_as_char = paste( "0", as.character(value), sep="" )
  }
  
  return(value_as_char)  
}

#using year-month-day-hour we can find the files where the data sought might be.
#this function creates
get_search_pattern_from_datetime <- function(date) {
  #file_name is arrival date, date must be shifted backwards 30 minutes.
  date = date + minutes(30)
  yr = as.character( year(date) )
  mo = char_and_add_zero( month(date) )
  dy = char_and_add_zero( day(date) )
  ho = char_and_add_zero( hour(date) )
  file_pattern = paste("c", yr, mo, dy, ho, sep="_")
  return(file_pattern)
}

get_filename_from_datetime <- function(date, or) {
  #file_name is arrival date, date must be shifted backwards 30 minutes.
  date = date + minutes(30)
  yr = as.character( year(date) )
  mo = char_and_add_zero( month(date) )
  dy = char_and_add_zero( day(date) )
  ho = char_and_add_zero( hour(date) )
  mi = char_and_add_zero( minute(date) )
  se = char_and_add_zero( second(date) )
  file_time = paste(ho,mi,se, sep = "")
  file_date = paste("c", yr, mo, dy, file_time, or, "cal.mat", sep="_")
  return(file_date)
}

date_difference_as_samples <- function(date_post, date_prev) {
  seconds_diff = as.numeric(date_post - date_prev, units = "secs")
  samples_diff = duration_to_samples(seconds_diff)
  return(samples_diff)
}

retrieve_file_list_by_orientation <- function(file_path, orientation) {
  search_path = dir( file_path )
  is_oriented = str_detect(search_path, orientation)
  is_filtered = str_detect(search_path, "fil")
  return( search_path[is_oriented & is_filtered] )
}

load_SR_data <- function(data_path, file) {
  data <- unlist( readMat( paste( data_path, file, sep="/" ) ), use.names = FALSE )
  return( data )
}

mean_na_handling <- function(values, threshold){
  #calculates the mean of the array given as values if the % of na is under threshold
  na_ind = is.na( values )
  na_ratio = sum( na_ind ) / length(values)
  
  if ( na_ratio <= threshold ) {
    mean_value = mean(values[ !na_ind ])
  } else {
    mean_value = NA
  }
  
  return(mean_value)
}

ma_with_na_handling <- function(values, threshold, w_before, w_after) { 
  #moving average using mean_na_handling; interpolates linearly any na elements
  #left after the filtering.
  averaged = slider::slide_dbl(values, mean, .before = w_before, .after = w_after, 
                               .complete = FALSE, .size = length( values ) )
  na_ind = which( is.na(averaged) )
  
  for (val in na_ind) {
    start_val = val - w_before;
    end_val = val + w_after;
    
    if ( end_val > length(values) ) {
      end_val = length(values)
    }
    
    if ( start_val >= 1 ) {
      averaged[val] = mean_na_handling( values[ start_val:end_val ], threshold )
    }
  }
  
  averaged = na.approx(averaged, na.rm = FALSE)
  return(averaged)
}

######### FUNCTIONS DEPENDENT ON GLOBAL VARIABLES ####

duration_to_samples <- function(duration) {
  return(duration * G_sample_frequency)
}

do_segment_exists <- function(date_start, orientation, segment_duration){
  full_path = paste(G_data_file_path,year(date_start),char_and_add_zero(month(date_start)),sep = "/")
  file_list = retrieve_file_list_by_orientation(full_path, orientation)
  file_pattern = get_search_pattern_from_datetime(date_start)
  I_candidate_files = which( str_detect(file_list, file_pattern) )
  #add the immediate previous segment to ensure all data from the queried hour
  #is selected
  if (length(I_candidate_files) == 0) {
    return(FALSE)
  }
  
  if (head(I_candidate_files, 1) != 1 ) { 
    I_candidate_files = c( head(I_candidate_files, 1) -1, I_candidate_files )
  }
  
  candidate_dates = lapply(file_list[I_candidate_files], get_datetime_from_file)
  I_chosen_from_candidates = max( which(candidate_dates < date_start) )
  I_chosen_file = I_candidate_files[I_chosen_from_candidates]
  
  sample_start = date_difference_as_samples(date_start,
                                            candidate_dates[[I_chosen_from_candidates]]
  ) + 1 #start at the beginning of the second.
  
  sample_end = date_difference_as_samples(date_start + seconds(segment_duration),
                                          candidate_dates[[I_chosen_from_candidates]]
  )
  
  if ( sample_end > G_sample_frequency*60*30 ) {
    if ( I_chosen_file == length(file_list) ) {
      return (FALSE)
    }      
    #concatenate the whole next file
  }
  
  return(TRUE)
}

is_data_interval_complete <- function(date_start, date_end, orientation, segment_duration, time_diff) {
  date_list = seq(date_start, date_end, by = time_diff)
  v_try_segment <- Vectorize(do_segment_exists, c("date_start") )
  
  segment_notif = v_try_segment(date_list, orientation, segment_duration)
  
  return(segment_notif)
}

load_segment_by_start_datetime <- function(date_start, orientation, segment_duration) {
  full_path = paste(G_data_file_path,year(date_start),char_and_add_zero(month(date_start)),sep = "/")
  file_list = retrieve_file_list_by_orientation(full_path, orientation)
  file_pattern = get_search_pattern_from_datetime(date_start)
  I_candidate_files = which( str_detect(file_list, file_pattern) )
  #add the immediate previous segment to ensure all data from the queried hour
  #is selected
  if (length(I_candidate_files) == 0) {
    warning(paste("no files found for ", date_start)) 
    return(NA)
  }
  
  if (head(I_candidate_files, 1) != 1 ) { 
    I_candidate_files = c( head(I_candidate_files, 1) -1, I_candidate_files )
  }
  
  candidate_dates = lapply(file_list[I_candidate_files], get_datetime_from_file)
  I_chosen_from_candidates = max( which(candidate_dates < date_start) )
  I_chosen_file = I_candidate_files[I_chosen_from_candidates]
  raw_data = load_SR_data( full_path, file_list[I_chosen_file] )
  
  if (G_debug) {
    print( paste("just read the data for ", date_start) )
  }
  
  sample_start = date_difference_as_samples(date_start,
                                            candidate_dates[[I_chosen_from_candidates]]
  ) + 1 #start at the beginning of the second.
  
  sample_end = date_difference_as_samples(date_start + seconds(segment_duration),
                                          candidate_dates[[I_chosen_from_candidates]]
  )
  
  if ( sample_end > length(raw_data) ) {
    if ( I_chosen_file == length(file_list) ) {
      warning(paste("Reached the end of file list for ", date_start)) 
      return (NA)
    }      
    #concatenate the whole next file
    if (G_debug) {
      print( paste("Gotta pick more for ", date_start) )
    }
    
    raw_data = c(raw_data, load_SR_data( full_path, file_list[I_chosen_file + 1] ))
  }
  
  segment_out = data = raw_data[ seq(sample_start, sample_end) ]
  
  return(segment_out)
}

load_segments_as_tsibble <- function(date_start, date_end, orientation, segment_duration, time_diff) {
  date_list = seq(date_start, date_end, by = time_diff)
  v_load_segment <- Vectorize(load_segment_by_start_datetime, c("date_start") )
  
  segment_matrix = v_load_segment(date_list, orientation, segment_duration)
  
  
  if (any(is.na(segment_matrix))) {
    #error propagation from load_segment_by_start_datetime
    warning("one or more files are missing!")
    if( is.matrix(segment_matrix) ) {
      return( rep( NA, ncol(segment_matrix) ) )
    }
    else{
      return( rep( NA, length(segment_matrix) ) )
    }
  }
  
  
  matrix_height = dim(segment_matrix)[2]  
  segment_data = vector("list", matrix_height)
  sg_mean = vector("numeric", matrix_height)
  sg_sd = vector("numeric", matrix_height)
  sg_skew = vector("numeric", matrix_height)
  sg_kurt = vector("numeric", matrix_height)
  
  for ( val in 1:matrix_height ) {
    segment_data[[val]] = segment_matrix[,val]
    sg_mean[val] = mean(segment_data[[val]])
    sg_sd[val] = sd(segment_data[[val]])
    sg_skew[val] = timeDate::skewness(segment_data[[val]])
    sg_kurt[val] = timeDate::kurtosis(segment_data[[val]], method = "moment")
    
  }
  
  data_frame <- tsibble(
    timestamp = as_datetime(date_list),
    data = segment_data,
    index = timestamp,
    mean = sg_mean,
    sd = sg_sd,
    skew = sg_skew,
    kurt = sg_kurt
  )
  
  return(data_frame)
}

perform_iterative_two_sample_test <- function(ts_segments, I_target) {
  current_date = as.character( ts_segments$timestamp[I_target] )
  target_segment = ts_segments$data[[ I_target ]]
  test_values = ts_segments %>%
    pmap_dbl( ~ks.test(..2, as.numeric(target_segment) )$statistic)
  
  ts_results = tsibble(
    timestamp = ts_segments$timestamp,
    ref_segment := current_date,
    KS_distance = test_values,
    index = timestamp,
    key = ref_segment
  )
  
  return(ts_results)
}

chain_iterative_two_sample_test <- function(ts_segments, I_targets) {
  ts_results = perform_iterative_two_sample_test(ts_segments, I_targets[1])
  
  if(length(I_targets) > 1){
    for (iter in 2:length(I_targets) ) {
      tmp_results = perform_iterative_two_sample_test(ts_segments, I_targets[iter])
      ts_results <- ts_results %>% rows_insert (tmp_results, by = c("timestamp","ref_segment") )
    }
  }
  
  return(ts_results)
}

downsample_array <- function(data_array_list, downsampling_factor) {
  data_array = unlist(data_array_list)
  downsampled_array = data_array[seq(1,length(data_array), by = downsampling_factor)]
  return(downsampled_array)
}

ma_smooth_hourly <- function(values, origin_pos) {
  averaged = ma_with_na_handling(values, na_ratio_threshold, 4, 5)
  averaged = ma_with_na_handling(values, na_ratio_threshold, 5, 4)
  averaged[origin_pos] = 0
  return(averaged)
}

fft_values <- function(signal, inverse) {
  signal_size = length(signal)
  signal_fft = fft(signal, inverse)
  signal_fft = abs(signal_fft/signal_size)
  signal_fft = signal_fft[1:(signal_size/2+1)]
  signal_fft[2:(length(signal_fft)-1)] = 2*signal_fft[2:(length(signal_fft)-1)]
  f = G_sample_frequency*seq(from = 0, to = signal_size/2)/signal_size
  signal_freq = tibble(
    freq_values = f,
    spectra_values = signal_fft
  ) 
  return(signal_freq)
}

averaged_fft <- function(t_segments, inverse){
  for (i_segments in 1:nrow(t_segments)) {
    tmp_signal = fft_values(t_segments$data[[i_segments]], inverse)
    if (i_segments == 1) {
      averaged_signal = tmp_signal
    }
    else {
      averaged_signal$spectra_values = averaged_signal$spectra_values + tmp_signal$spectra_values
    }
  }
  averaged_signal$spectra_values = averaged_signal$spectra_values/i_segments
  return(averaged_signal)
}
