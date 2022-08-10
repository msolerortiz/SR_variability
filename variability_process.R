############ Freq Graphs ####
date_start = ymd_hms("2018-04-20 10:00:00")
date_end = ymd_hms("2018-04-20 11:00:00")
forward_step = 7
duration = 7

t_KS_distances = tibble (
  datepoint = ymd_hms(),
  KS_distance = numeric(),
  c_year = numeric(),
  c_month = numeric(),
  c_day = numeric(),
  c_hour = numeric(),
  c_minute = numeric(),
  c_second = numeric()
)

ts_segment = load_segments_as_tsibble( date_start, date_end, "NS", duration, forward_step )
if ( !any(is.na(ts_segment) ) ) {
  #ts_segment = ts_segment %>%
  #mutate( data = map( .x = data, ~downsample_array(.x, downsampling_factor = 50) ) )
  
  target_segment = ceiling(nrow(ts_segment)/2)
  
  ts_tmp = chain_iterative_two_sample_test( ts_segment, target_segment )
  
  target_date = ts_tmp$timestamp[target_segment]
  
  total_datapoints = nrow(ts_segment)
  
  t_KS_distances = t_KS_distances %>%
    add_row( datepoint = ts_tmp$timestamp,
             KS_distance = ts_tmp$KS_distance,
             c_year = rep( year( target_date ), total_datapoints ), 
             c_month = rep( month( target_date ), total_datapoints ), 
             c_day = rep( day( target_date ), total_datapoints ), 
             c_hour = rep( hour( target_date ), total_datapoints ),
             c_minute = rep( minute( target_date ), total_datapoints ),
             c_second = rep( second( target_date ), total_datapoints )
    )
  
}

t_KS_distances %>%
  filter(KS_distance != -1) %>%
  #group_by(duration) %>%
  #mutate( KS_distance = ma_smooth_hourly(KS_distance,target_segment) ) %>%
  ggplot( aes( x= datepoint, y = KS_distance ) ) +
  geom_line(size = 0.8) +
  scale_x_datetime(name = "Datetime",
  )+
  scale_y_continuous(name = "KS distance") +
  geom_hline(yintercept = 0.03, color = "#D39200", size = 0.8) +
  geom_hline(yintercept = 0.11, color = "#00B9E3", size = 0.8)

fft_segments_low = ts_segment[which(t_KS_distances$KS_distance < 0.03),]
fft_segments_high = ts_segment[which(t_KS_distances$KS_distance > 0.11),]

avg_fft_low = averaged_fft(fft_segments_low, FALSE) %>%
  mutate(freq_values = round(freq_values, digits = 4) )
avg_fft_high = averaged_fft(fft_segments_high, FALSE) %>%
  mutate(freq_values = round(freq_values, digits = 4) )

num_samples_std = ( nrow(fft_segments_low) + nrow(fft_segments_high) )/2

avg_fft_std = averaged_fft(ts_segment[1:num_samples_std,], FALSE) %>%
  mutate(freq_values = round(freq_values, digits = 4) )

avg_fft = avg_fft_low %>%
  full_join(avg_fft_high, by ="freq_values", suffix = c(".low",".high") ) %>%
  full_join(avg_fft_std, by = "freq_values" ) %>%
  pivot_longer(c(spectra_values.low, spectra_values.high, spectra_values), names_to = "key", values_to = "value")



avg_fft %>%
  filter(freq_values > 4, freq_values < 50) %>%
  ggplot(aes(x = freq_values, y = value, color = as.factor(key) ) ) +
  geom_line( size = 1.2 )+
  scale_colour_manual(
    values = c("#F8766D", "#00B9E3", "#D39200", "#00C19F", "#00BA38"),
    name = "Segments chosen",
    labels = c("First segments", "Divergent segments", "Similar segments" )
  ) + 
  scale_x_continuous(name = "Frequency [Hz]",
  )+
  scale_y_continuous(name = expression(paste("Intensity [pT/",sqrt(Hz),"]") ) )

############ DATA PROCESS ####
date_start = ymd_hms("2019-07-01 00:00:00")
date_end = ymd_hms("2019-08-01 00:00:00")
forward_step = 60*2
duration = 60*10
time_interval = 60*60 #in seconds

s_dates = seq(from = date_start, to = date_end, by = as.difftime(time_interval, units = "secs") )

t_KS_distances = tibble (
  n_datapoint = numeric(),
  KS_distance = numeric(),
  c_year = numeric(),
  c_month = numeric(),
  c_day = numeric(),
  c_hour = numeric()
)

for (d_iter in 2:length(s_dates) ) {
  ts_segment = load_segments_as_tsibble( s_dates[d_iter - 1], s_dates[ d_iter ], "EW", duration, forward_step )
  # if load_segments_as_tsibble misses any amount of data, it will return as much NAs
  # as rows the tsibble would have if it had data.
  
  if ( any(is.na(ts_segment) ) ){
    insert_size = length(ts_segment)
    distances_array = rep( -1, insert_size )
    target_segment = ceiling(insert_size/2) #TARGET SEGMENT IS ALWAYS THE MIDDLE OF THE INTERVAL   
  }
  
  else {
    ts_segment = ts_segment %>%
      mutate( data = map( .x = data, ~downsample_array(.x, downsampling_factor = 100) ) )
    insert_size = nrow(ts_segment)
    target_segment = ceiling(insert_size/2)          
    ts_tmp = chain_iterative_two_sample_test( ts_segment, target_segment )
    
    distances_array = ts_tmp$KS_distance
  }
  
  index_values = seq(from = 1, to = insert_size, by = 1) -  target_segment
  
  t_KS_distances = t_KS_distances %>% add_row(n_datapoint = index_values, 
                                              KS_distance = distances_array,
                                              c_year = rep( year( s_dates[d_iter - 1] ), insert_size ), 
                                              c_month = rep( month( s_dates[d_iter - 1] ), insert_size ), 
                                              c_day = rep( day( s_dates[d_iter - 1] ), insert_size ), 
                                              c_hour = rep( hour( s_dates[d_iter - 1] ), insert_size ) 
  )
}


t_KS_distances_all %>%
  filter(KS_distance != -1) %>%
  group_by(c_year, c_month, c_hour) %>%
  summarize(
    mean_value = mean(KS_distance),
    std_value = sd(KS_distance),
    .groups = "keep"
  ) %>%
  pivot_longer(c(mean_value, std_value), names_to = "key", values_to = "value") %>%
  unite(key_year, c("key", "c_year")) %>%
  ggplot( aes( x= c_hour, y = c(value), colour = as.factor(key_year) ) ) +
  facet_wrap(~ month(c_month, abbr = FALSE, label = TRUE), nrow = 3 ) +
  geom_line()+
  scale_colour_manual(
    values = c("springgreen", "darkolivegreen4", "lightskyblue", "lightskyblue4"),
    name = "Statistic and year",
    labels = c("Mean 2018", "Mean 2019", "Std. deviation 2018", "Std. deviation 2019")
  ) + 
  scale_x_continuous(name = "Time [UTC]",
                     breaks = c(0,4,8,12,16,20)
  )+
  scale_y_continuous(name = "KS distance")+
  theme(legend.position = "top", 
        legend.box = "horizontal")



############ PARAMETER TEST - Duration####
date = ymd_hms("2018-12-27 00:00:00")
forward_step = 60*2
duration = c(60, 60*1.5, 60*2, 60*5, 60*10)
time_interval = days(1)
t_KS_distances = tibble (
  datepoint = ymd_hms(),
  KS_distance = numeric(),
  duration = numeric(),
  c_year = numeric(),
  c_month = numeric(),
  c_day = numeric(),
  c_hour = numeric(),
  c_minute = numeric(),
  c_second = numeric()
)

if ( !( any( is.na( is_data_interval_complete( date, date + time_interval, "EW", duration[1], forward_step ) ) ) ) ) {
  
  for (iter in 1:length(duration)) {
    ts_segment = load_segments_as_tsibble( date, date + time_interval, "EW", duration[iter], forward_step )
    ts_segment = ts_segment %>%
      mutate( data = map( .x = data, ~downsample_array(.x, downsampling_factor = 20*iter) ) )
    
    target_segment = ceiling(nrow(ts_segment)/2)
    
    ts_tmp = chain_iterative_two_sample_test( ts_segment, target_segment )
    
    target_date = ts_tmp$timestamp[target_segment]
    
    total_datapoints = nrow(ts_segment)
    
    t_KS_distances = t_KS_distances %>%
      add_row( datepoint = ts_tmp$timestamp,
               KS_distance = ts_tmp$KS_distance,
               duration = rep( duration[iter], total_datapoints ),
               c_year = rep( year( target_date ), total_datapoints ), 
               c_month = rep( month( target_date ), total_datapoints ), 
               c_day = rep( day( target_date ), total_datapoints ), 
               c_hour = rep( hour( target_date ), total_datapoints ),
               c_minute = rep( minute( target_date ), total_datapoints ),
               c_second = rep( second( target_date ), total_datapoints )
      )
    
  }
  
}

t_KS_distances %>%
  filter(KS_distance != -1) %>%
  group_by(duration) %>%
  mutate( KS_distance = ma_smooth_hourly(KS_distance,target_segment) ) %>%
  ggplot( aes( x= datepoint, y = KS_distance, colour = as.factor(duration) ) ) +
  geom_line(size = 0.6) +
  scale_colour_manual(
    values = c("#F8766D", "#00B9E3", "#D39200", "#800080", "#00BA38"),
    name = "Segment duration [s]",
  ) + 
  scale_x_datetime(name = "Datetime [UTC]",
                   date_breaks = "2 hours",
                   date_labels = "%H:%M"
  )+
  scale_y_continuous(name = "KS distance")+
  theme(legend.position = "top",
        legend.box = "horizontal")


############ PARAMETER TEST - Forward Step ####
date = ymd_hms("2018-06-24 15:00:00")
forward_step = c(60/30, 60/10, 60/5, 60)
duration = 60*5
time_interval = hours(1)

t_KS_distances = tibble (
  datepoint = ymd_hms(),
  KS_distance = numeric(),
  forward_step = numeric(),
  c_year = numeric(),
  c_month = numeric(),
  c_day = numeric(),
  c_hour = numeric(),
  c_minute = numeric(),
  c_second = numeric()
)

if ( !( any( is.na( is_data_interval_complete( date, date + time_interval, "EW", duration, forward_step[1] ) ) ) ) ) {
  
  for (iter in 1:length(forward_step)) {
    ts_segment = load_segments_as_tsibble( date, date + time_interval, "EW", duration, forward_step[iter] )
    ts_segment = ts_segment %>%
      mutate( data = map( .x = data, ~downsample_array(.x, downsampling_factor = 50) ) )
    
    target_segment = ceiling(nrow(ts_segment)/2)
    
    ts_tmp = chain_iterative_two_sample_test( ts_segment, target_segment )
    
    target_date = ts_tmp$timestamp[target_segment]
    
    total_datapoints = nrow(ts_segment)
    
    t_KS_distances = t_KS_distances %>%
      add_row( datepoint = ts_tmp$timestamp,
               KS_distance = ts_tmp$KS_distance,
               forward_step = rep( forward_step[iter], total_datapoints ),
               c_year = rep( year( target_date ), total_datapoints ), 
               c_month = rep( month( target_date ), total_datapoints ), 
               c_day = rep( day( target_date ), total_datapoints ), 
               c_hour = rep( hour( target_date ), total_datapoints ),
               c_minute = rep( minute( target_date ), total_datapoints ),
               c_second = rep( second( target_date ), total_datapoints )
      )
    
  }
  
}

t_KS_distances %>%
  filter(KS_distance != -1) %>%
  group_by(forward_step) %>%
  #mutate( KS_distance = ma_smooth_hourly(KS_distance,target_segment) ) %>%
  ggplot( aes( x= datepoint, y = KS_distance, colour = as.factor(forward_step) ) ) +
  geom_line() +
  scale_colour_manual(
    values = c("#F8766D", "#00B9E3", "#D39200", "#00C19F", "#00BA38"),
    name = "Segment duration",
  ) + 
  scale_x_datetime(name = "Datetime",
  )+
  scale_y_continuous(name = "KS distance")

t_KS_distances %>%
  filter(KS_distance != -1, forward_step == 2) %>%
  mutate( smooth_KS_distance = ma_smooth_hourly(KS_distance,nrow(.)/2) ) %>%
  pivot_longer(c(KS_distance, smooth_KS_distance), names_to = "key", values_to ="value") %>%
  ggplot( aes( x= datepoint, y = value, colour = as.factor(key) ) ) +
  geom_line() +
  scale_colour_manual(
    values = c("#F8766D", "#00B9E3", "#D39200", "#00C19F", "#00BA38"),
    name = NULL,
    labels = c("Raw", "Smoothed")
  ) + 
  scale_x_datetime(name = "Time [UTC]",
                   date_labels = "%H:%M"
  )+
  scale_y_continuous(name = "KS distance")+
  theme(legend.position = "top",
        legend.box = "horizontal")


############ DIFFERENT SEGMENTS, SAME DAY ####
date = ymd_hms("2018-02-04 00:00:00")
forward_step = 60*2
duration = 60*10
time_interval = days(1)

t_KS_distances = tibble (
  datepoint = ymd_hms(),
  KS_distance = numeric(),
  position = numeric(),
  c_year = numeric(),
  c_month = numeric(),
  c_day = numeric(),
  c_hour = numeric(),
  c_minute = numeric()
)

ts_segment = load_segments_as_tsibble( date, date + time_interval, "EW", duration, forward_step )
number_of_target_segments = 5

if ( !( any( is.na(ts_segment) ) ) ) {
  ts_segment = ts_segment %>%
    mutate( data = map( .x = data, ~downsample_array(.x, downsampling_factor = 100) ) )
  
  for (iter in 1:number_of_target_segments) {
    #target_segment = round( nrow(ts_segment)/2 )
    #target_segment = iter*round( nrow(ts_segment)/number_of_target_segments )
    target_segment = (iter-1)*round( nrow(ts_segment)/number_of_target_segments ) + 1
    
    target_date = ts_segment$timestamp[target_segment]
    
    ts_tmp = chain_iterative_two_sample_test( ts_segment, target_segment )
    
    total_datapoints = nrow(ts_segment)
    
    t_KS_distances = t_KS_distances %>%
      add_row( datepoint = ts_tmp$timestamp,
               KS_distance = ts_tmp$KS_distance,
               position = rep( target_segment, total_datapoints ),
               c_year = rep( year( target_date ), total_datapoints ), 
               c_month = rep( month( target_date ), total_datapoints ), 
               c_day = rep( day( target_date ), total_datapoints ), 
               c_hour = rep( hour( target_date ), total_datapoints ),
               c_minute = rep( minute( target_date ), total_datapoints )
      )
    
  }
  
  table_dates = t_KS_distances %>% 
    select(c_year,c_month,c_day,c_hour,c_minute) %>% 
    distinct() 
  
  target_segment_dates = ymd_hm( paste(
    paste(table_dates$c_year, table_dates$c_month, table_dates$c_day, sep = "-"),
    paste(table_dates$c_hour, table_dates$c_minute, sep=":"),
    sep = " " )
  )
  
  target_segment_dates = rep(target_segment_dates, each = total_datapoints)
  
}
#################### plot for many dates
t_KS_distances %>%
  filter(KS_distance != -1) %>%
  group_by(c_year, c_month, c_day, c_hour, c_minute) %>%
  mutate( KS_distance = ma_smooth_hourly(KS_distance,position) ) %>%
  ggplot( aes( x = datepoint, 
               y = KS_distance, 
               colour = as.factor(format(target_segment_dates, "%H:%M") ) ) ) +
  geom_line(size = 0.6) +
  scale_colour_manual(
    values = c("#F8766D", "#00B9E3", "#D39200", "#800080", "#00BA38"),
    name = "Target segment time [UTC]",
  ) + 
  scale_x_datetime(name = "Datetime [UTC]",
                   date_breaks = "2 hours",
                   date_labels = "%H:%M"
  )+
  scale_y_continuous(name = "KS distance")+
  theme(legend.position = "top",
        legend.box = "horizontal")



#################### plot to show raw vs smoothed
t_KS_distances %>%
  filter(KS_distance != -1) %>%
  group_by(c_year, c_month, c_day, c_hour, c_minute) %>%
  mutate( smooth_KS_distance = ma_smooth_hourly(KS_distance,position) ) %>%
  pivot_longer(c(KS_distance, smooth_KS_distance), names_to = "key", values_to ="value") %>%
  ggplot( aes( x = datepoint, 
               y = value, 
               colour = as.factor(key) ) ) +
  geom_line() +
  scale_colour_manual(
    values = c("#F8766D", "#00B9E3", "#D39200", "#00C19F", "#00BA38"),
    name = "Target segment",
  ) + 
  scale_x_datetime(name = "Datetime",
  )+
  scale_y_continuous(name = "KS distance")