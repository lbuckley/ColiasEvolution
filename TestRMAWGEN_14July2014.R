#ComprehensiveTemperatureGenerator

station = c("T0001", "T0010", "T0099"); Tx_all; Tn_all; 
    mean_climate_Tn = NULL; mean_climate_Tx = NULL; Tx_spline = NULL; 
    Tn_spline = NULL; year_max = 1990; year_min = 1961; leap = TRUE; 
    nmonth = 12; verbose = TRUE; p = 1; type = "none"; lag.max = NULL; 
    ic = "AIC"; activateVARselect = FALSE; year_max_sim = year_max; 
    year_min_sim = year_min; mean_climate_Tn_sim = NULL; mean_climate_Tx_sim = NULL; 
    Tn_spline_sim = NULL; Tx_spline_sim = NULL; onlygeneration = FALSE; 
    varmodel = NULL; normalize = TRUE; type_quantile = 3; sample = NULL; 
    extremes = TRUE; option = 2; yearly = FALSE; yearly_sim = yearly;
    n_GPCA_iteration = 0; n_GPCA_iteration_residuals = n_GPCA_iteration; 
    exogen = NULL; exogen_sim = exogen; is_exogen_gaussian = FALSE; 
    exogen_all = NULL; exogen_all_col = station; nscenario = 1; 
    seed = NULL; noise = NULL

##TEST
station=vstation; Tx_all=TEMPERATURE_MAX; Tn_all=TEMPERATURE_MIN; year_min=year_min; year_max=year_max
 p=p; n_GPCA_iteration=n_GPCA_iter; n_GPCA_iteration_residuals=n_GPCA_iteration_residuals
 sample="monthly"; year_min_sim=year_min_sim; year_max_sim=year_max_sim


    useVAR = TRUE
    if (option == 2) 
        normalize = TRUE
    origin <- paste(year_min, "1", "1", sep = "/")
    origin_sim <- paste(year_min_sim, "1", "1", sep = "/")
    param <- setComprehensiveTemperatureGeneratorParameters(station = station, 
        Tx_all = Tx_all, Tn_all = Tn_all, mean_climate_Tn = mean_climate_Tn, 
        mean_climate_Tx = mean_climate_Tx, Tx_spline = Tx_spline, 
        Tn_spline = Tn_spline, year_max = year_max, year_min = year_min, 
        leap = leap, nmonth = nmonth, verbose = verbose, cpf = NULL, 
        normalize = normalize, sample = sample, option = option, 
        yearly = yearly)

#---------------
#setComprehensiveTemperatureGeneratorParameters
#function (station, Tx_all, Tn_all, 
mean_climate_Tn = NULL; mean_climate_Tx = NULL; 
Tx_spline = NULL; Tn_spline = NULL; year_max = 1990; year_min = 1961; 
    leap = TRUE; nmonth = 12; verbose = FALSE; cpf = NULL; normalize = TRUE; 
    sample = NULL; option = 2; yearly = FALSE 
  
    station <- station[(station %in% names(Tx_all)) & (station %in% 
        names(Tn_all))]
    origin <- paste(year_min, "1", "1", sep = "-")
    Tn_mes <- as.data.frame(extractyears(Tn_all, year_min = year_min, 
        year_max = year_max, station = station))
#-----------------
#Fails to extract data
as.data.frame(extractyears(Tn_all, year_min = 1953, year_max = 1954, station = c("C1","D1")) )
names(Tn_all)= c("year","month","day","A1","B1","C1","D1")


#-----------------
    Tx_mes <- as.data.frame(extractyears(Tx_all, year_min = year_min, 
        year_max = year_max, station = station))
    names(Tn_mes) <- station
    names(Tx_mes) <- station
    Tm_mes <- (Tn_mes + Tx_mes)/2
    DeltaT_mes <- Tx_mes - Tn_mes
    if (!is.monthly.climate(mean_climate_Tx, nstation = length(station), 
        nmonth = nmonth, verbose = verbose)) 
        mean_climate_Tx <- getMonthlyMean(as.data.frame(Tx_mes), 
            year_min = year_min, year_max = year_max, station = station, 
            no_date = TRUE, origin = origin, yearly = yearly)
#---------------------
#getMonthlyMean

function (data, year_min = 1961, year_max = 1990, station = names(data), 
    no_date = FALSE, origin = "1961-1-1", yearly = FALSE) 
{
    if (is.null(station)) {
        nstation = ncol(data)
    }
    else {
        nstation = length(station)
    }
    if (no_date) {
        dates <- findDate(1:nrow(data), origin = origin, data.frame = TRUE, 
            decimal = FALSE, character = FALSE)
        data <- cbind(dates, data)
        datalocked <- data
    }
    if (yearly) {
        out <- list()
        for (year in year_min:year_max) {
            out[[as.character(year)]] <- getMonthlyMean(data, 
                year_min = year, year_max = year, station = station, 
                no_date = FALSE, origin = NULL, yearly = FALSE)
        }
    }
    else {
        nmonth = 12
        out <- as.data.frame(array(NA, c(nmonth, nstation)))
        names(out) <- station
        for (r in 1:nmonth) {
            im <- (data$year <= year_max) & (data$year >= year_min) & 
                (data$month == r)
            for (c in 1:nstation) {
                var <- data[im, station[c]]
                out[r, c] <- mean(var[!is.na(var)])
            }
        }
    }
    return(out)
}

#----------------------
    if (!is.monthly.climate(mean_climate_Tn, nstation = length(station), 
        nmonth = nmonth, verbose = verbose)) 
        mean_climate_Tn <- getMonthlyMean(as.data.frame(Tn_mes), 
            year_min = year_min, year_max = year_max, station = station, 
            no_date = TRUE, origin = origin, yearly = yearly)
    monthly_mean_Tx <- mean_climate_Tx
    monthly_mean_Tn <- mean_climate_Tn
    nyear <- year_max - year_min + 1
    if (is.null(Tx_spline)) {
        Tx_spline <- as.data.frame(splineInterpolateMonthlytoDailyforSeveralYears(val = mean_climate_Tx, 
            start_year = year_min, nyear = nyear, leap = leap, 
            yearly = yearly))
        if (yearly) {
            names(Tx_spline) <- colnames(mean_climate_Tx[[1]])
        }
        else {
            names(Tx_spline) <- colnames(mean_climate_Tx)
        }
    }
    if (is.null(Tn_spline)) {
        Tn_spline <- as.data.frame(splineInterpolateMonthlytoDailyforSeveralYears(val = mean_climate_Tn, 
            start_year = year_min, nyear = nyear, leap = leap, 
            yearly = yearly))
        if (yearly) {
            names(Tn_spline) <- colnames(mean_climate_Tn[[1]])
        }
        else {
            names(Tn_spline) <- colnames(mean_climate_Tn)
        }
    }
    SplineAdvTm <- (Tx_spline + Tn_spline)/2
    SplineAdvDeltaT <- (Tx_spline - Tn_spline)
    nstation = length(station)
    Tn_mes_res <- Tn_mes - Tn_spline
    Tx_mes_res <- Tx_mes - Tx_spline
    Tm_mes_res <- Tm_mes - SplineAdvTm
    DeltaT_mes_res <- DeltaT_mes
    stdTn <- apply(X = Tn_mes_res, MARGIN = 2, FUN = sd, na.rm = TRUE)
    stdTx <- apply(X = Tx_mes_res, MARGIN = 2, FUN = sd, na.rm = TRUE)
    stdTm <- apply(X = Tm_mes_res, MARGIN = 2, FUN = sd, na.rm = TRUE)
    for (s in 1:nstation) {
        Tn_mes_res[, s] <- (Tn_mes[, s] - Tn_spline[, s])/stdTn[s]
        Tx_mes_res[, s] <- (Tx_mes[, s] - Tx_spline[, s])/stdTx[s]
        Tm_mes_res[, s] <- (Tm_mes[, s] - SplineAdvTm[, s])/stdTm[s]
    }
    if (option == 1) {
        data_original <- cbind(Tx_mes_res, Tn_mes_res)
    }
    else if (option == 2) {
        data_original <- cbind(Tm_mes_res, DeltaT_mes_res)
        normalize = TRUE
    }
    else if (option == 3) {
        data_original <- Tm_mes_res
    }
    else {
        data_original <- NULL
    }
    if (normalize) {
        data_for_var <- normalizeGaussian_severalstations(x = data_original, 
            data = data_original, sample = sample, cpf = cpf, 
            origin_x = origin, origin_data = origin)
    }
    else {
        data_for_var <- data_original
    }
    out <- list(data_for_var, data_original, Tn_mes_res, Tx_mes_res, 
        Tm_mes_res, DeltaT_mes_res, stdTn, stdTx, stdTm, Tn_spline, 
        Tx_spline, SplineAdvTm, SplineAdvDeltaT, Tn_mes, Tx_mes, 
        Tm_mes, DeltaT_mes, monthly_mean_Tx, monthly_mean_Tn)
    names(out) <- c("data_for_var", "data_original", "Tn_mes_res", 
        "Tx_mes_res", "Tm_mes_res", "DeltaT_mes_res", "stdTn", 
        "stdTx", "stdTm", "Tn_spline", "Tx_spline", "SplineAdvTm", 
        "SplineAdvDeltaT", "Tn_mes", "Tx_mes", "Tm_mes", "DeltaT_mes", 
        "monthly_mean_Tx", "monthly_mean_Tn")



#---------------    
if (!onlygeneration) {
        if (!is.null(exogen_all)) {
            exogen <- as.data.frame(extractyears(exogen_all, 
                year_min = year_min, year_max = year_max, station = exogen_all_col))
            is_exogen_gaussian = FALSE
            if (is.null(exogen_sim)) 
                exogen_sim <- exogen
        }
        if (!is.null(exogen) & (!is_exogen_gaussian)) {
            exogen0 <- exogen
            exogen <- normalizeGaussian_severalstations(x = exogen0, 
                data = exogen0, sample = sample, cpf = NULL, 
                origin_x = origin, origin_data = origin, extremes = extremes)
        }
        var <- getVARmodel(data = param[["data_for_var"]], suffix = c("_T1", 
            "_T2"), sep = "", p = p, type = type, lag.max = lag.max, 
            ic = ic, activateVARselect = activateVARselect, exogen = exogen, 
            n_GPCA_iteration_residuals = n_GPCA_iteration_residuals, 
            n_GPCA_iteration = n_GPCA_iteration, extremes = extremes)
        if (activateVARselect) 
            return(list(input = param, varselect = var))
    }
    else {
        var <- varmodel
    }
    if (is.null(Tx_spline)) 
        Tx_spline <- param[["Tx_spline"]]
    if (is.null(Tn_spline)) 
        Tn_spline <- param[["Tn_spline"]]
    nyear_sim <- year_max_sim - year_min_sim + 1
    if (is.null(Tx_spline_sim)) {
        if (is.null(mean_climate_Tx_sim)) 
            mean_climate_Tx_sim <- param[["monthly_mean_Tx"]]
        Tx_spline_sim <- as.data.frame(splineInterpolateMonthlytoDailyforSeveralYears(val = mean_climate_Tx_sim, 
            start_year = year_min_sim, nyear = nyear_sim, leap = leap, 
            yearly = yearly_sim))
        if (yearly_sim) {
            names(Tx_spline_sim) <- colnames(mean_climate_Tx_sim[[1]])
        }
        else {
            names(Tx_spline_sim) <- colnames(mean_climate_Tx_sim)
        }
    }
    if (is.null(Tn_spline_sim)) {
        if (is.null(mean_climate_Tn_sim)) 
            mean_climate_Tn_sim <- param[["monthly_mean_Tn"]]
        Tn_spline_sim <- as.data.frame(splineInterpolateMonthlytoDailyforSeveralYears(val = mean_climate_Tn_sim, 
            start_year = year_min_sim, nyear = nyear_sim, leap = leap, 
            yearly = yearly_sim))
        if (yearly_sim) {
            names(Tn_spline_sim) <- colnames(mean_climate_Tn_sim[[1]])
        }
        else {
            names(Tn_spline_sim) <- colnames(mean_climate_Tn_sim)
        }
    }
    if (!is.null(noise)) {
        if (noise == "residuals") 
            noise <- residuals(var)
    }
    SplineAdvTm_sim <- (Tx_spline_sim + Tn_spline_sim)/2
    SplineAdvDeltaT_sim <- (Tx_spline_sim - Tn_spline_sim)
    if (!is.null(exogen_sim) & (!is_exogen_gaussian)) {
        exogen0_sim <- exogen_sim
        exogen_sim <- normalizeGaussian_severalstations(x = exogen0_sim, 
            data = exogen0_sim, sample = sample, cpf = NULL, 
            origin_x = origin_sim, origin_data = origin_sim, 
            extremes = extremes)
    }
    if (!is.null(seed)) 
        set.seed(seed)
    original_data <- param[["data_original"]]
    if (option == 2) {
        ntall <- as.integer(ncol(original_data))
        ntn <- as.integer(ncol(original_data)/2) + 1
        if (nrow(Tx_spline_sim) < nrow(original_data)) {
            nfrac <- as.integer(nrow(original_data)/nrow(Tx_spline_sim)) + 
                1
            Tx_spline_sim2 <- mapply(rep, Tx_spline_sim, nfrac)
            Tn_spline_sim2 <- mapply(rep, Tn_spline_sim, nfrac)
        }
        else {
            Tx_spline_sim2 <- Tx_spline_sim
            Tn_spline_sim2 <- Tn_spline_sim
        }
        ntn_rows <- 1:nrow(original_data)
        original_data[, ntn:ntall] <- original_data[, ntn:ntall]/(Tx_spline[ntn_rows, 
            station] - Tn_spline[ntn_rows, station]) * (Tx_spline_sim2[ntn_rows, 
            station] - Tn_spline_sim2[ntn_rows, station])
    }
    results <- generateTemperatureTimeseries(std_tn = param[["stdTn"]], 
        std_tx = param[["stdTx"]], SplineTx = Tx_spline_sim, 
        SplineTn = Tn_spline_sim, SplineTm = SplineAdvTm_sim, 
        SplineDeltaT = SplineAdvDeltaT_sim, std_tm = param[["stdTm"]], 
        var = var, normalize = normalize, type = type_quantile, 
        sample = sample, option = option, original_data = original_data, 
        origin_x = origin_sim, origin_data = origin, exogen = exogen_sim, 
        extremes = extremes, noise = noise)
    if (nscenario > 1) {
        for (kk in 2:nscenario) {
            results_temp <- generateTemperatureTimeseries(std_tn = param[["stdTn"]], 
                std_tx = param[["stdTx"]], SplineTx = Tx_spline_sim, 
                SplineTn = Tn_spline_sim, SplineTm = SplineAdvTm_sim, 
                SplineDeltaT = SplineAdvDeltaT_sim, std_tm = param[["stdTm"]], 
                var = var, normalize = normalize, type = type_quantile, 
                sample = sample, option = option, original_data = param[["data_original"]], 
                origin_x = origin_sim, origin_data = origin, 
                exogen = exogen_sim, extremes = extremes)
            Tx_index <- sprintf("Tx_gen%05d", kk)
            Tn_index <- sprintf("Tn_gen%05d", kk)
            results[[Tx_index]] <- results_temp$Tx_gen
            results[[Tn_index]] <- results_temp$Tn_gen
        }
    }
    if (onlygeneration) {
        return(list(output = results))
    }
    else {
        return(list(input = param, var = var, output = results, 
            temporary = original_data))
    }
    return(list(input = param, var = var, output = results))
}
<environment: namespace:RMAWGEN>