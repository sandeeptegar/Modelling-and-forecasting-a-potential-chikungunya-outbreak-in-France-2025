# Writing function chikungunya()
# input arguements: latitude, longitude, t_start, t_end, seed_date
# latitude and longitides are in Float64 foramt while t_start and t_end are in date format "dd.mm.yyyy".
function chikungunya(latitude, longitude, t_start::String, t_end::String, seed_date::String)

    # Function to convert integer dates to date time format
    # example: x = 3.45 and datetime = DateTime(2013,7,1,12,30)
    function add_float2datetime(x,datetime)
        d = Int(floor(x))
        e = x-d
        h = Int(floor(e*24))
        e = e-h/24
        m = Int(floor(e*24*60))
        e = e-m/24/60
        s = Int(floor(e*24*60*60))
        e = e-s/24/60/60
        ms = Int(floor(e*24*60*60*1e3))
        return datetime+Day(d)+Hour(h)+Minute(m)+Second(s)+Millisecond(ms)
    end
    # Temperature dependent biting rate
    # G_Time(Temp) is the length of gonotrophic cycle.
    function bite(Temp)
        1 ./G_Time(Temp)
    end
    # Probability of infection from human to vector per bite
    # 1 Liu-Helmersson et al. 2014
    # Based on empirical data for various flaviviruses (West Nile virus, Murray Valley encephalitis virus,
    # and St. Louis encephalitis virus), Lambrechts et al 2011.
    #function human_vector(Temp)
    #    out = ifelse.(Temp .<= 26.1, 0.0729.*Temp .- 0.9037,1.0)
    #    if out <= 0
    #        out = 0
    #    end
    #    if Temp > 32.461
    #        out = 0
    #    end
    #    return(out)
    #end
    #2 Mordecai et al. 2017 (Probability of infection based on DENV-Ae. albopictus combination)
    function human_vector(Temp)
        if Temp > 36.82
            Temp = 36.82
        end
        out = 4.394e-4 .*Temp.*(Temp .- 3.62) .* (36.82.-Temp).^0.5
        if out <= 0
            out = 0
        end
        if out >= 1
            out = 1
        end
        if Temp < 5
            out = 0
        end
        return(out)
    end
    # The probability of transmission from vector to human per bite
    # for DENV
    #function vector_human(Temp)
     #   if Temp > 32.461
      #      Temp = 32.461
       # end
        #out = 0.7 .* 0.001044.*Temp.*(Temp .- 12.286) .* (32.461.-Temp).^0.5
        #if out <= 0
         #   out = 0
        #end
        #if out >= 1
         #   out = 1
        #end
        #if Temp < 5
         #   out = 0
        #end
        #return(out)
    #end
    # Modified according to our model for chikungunya (quadratic function)
    function vector_human(Temp)
        if Temp > 34.73
            Temp = 34.73
        end
        out = 1.0.*0.0037.*(Temp .- 8.93) .* (34.73.-Temp)
        if out <= 0
            out = 0
        end
        if out >= 1
            out = 1
        end
        if Temp < 8.93
            out = 0
        end
        return(out)
    end
    
    # Loads climates data
    filename = string(latitude, "_", longitude, ".csv")
    Guan_DP_Clim = CSV.read(filename, header=true, DataFrame)

    # read time index (use +1 as header is the first line in the csv)
    idx_1 = findfirst(==(Date(t_start, "d.m.yyyy")),Guan_DP_Clim.time)
    idx_2 = findfirst(==(Date(t_end, "d.m.yyyy")),Guan_DP_Clim.time)
    
    #Count the number of days between dates including end dates
    num_days = abs(Dates.value(Guan_DP_Clim.time[idx_2] - Guan_DP_Clim.time[idx_1]))

    # Reads in climate data and converts to splines (note DP doesn't do anything anymore)
    Guan_Temp = Guan_DP_Clim[:,4]
    Guan_Prec = Guan_DP_Clim[:,5]
    Guan_EV   = Guan_DP_Clim[:,6]
    Precip_spl = Spline1D(collect(range(1,num_days+1,length = num_days+1)), Guan_Prec[idx_1+1:idx_2+1])
    Temps_spl = Spline1D(collect(range(1,num_days+1,length = num_days+1)), Guan_Temp[idx_1+1:idx_2+1])
    EV_spl = Spline1D(collect(range(1,num_days+1,length = num_days+1)), Guan_EV[idx_1+1:idx_2+1]);
    # Extrinsic incubation period
    # for dengue
    #function EIP(Temp)
    #    out = exp(2.9 - 0.08 * Temp + 1/(2*4.9))
    #    if Temp < 12
    #        out = exp(2.9 - 0.08 * 12 + 1/(2*4.9))
    #    end
    #    return(out)
    #end
    # for chikungunya
    function EIP(Temp)
        out = 55.03*exp(-0.13*Temp/(1+(Temp/273)))
        if Temp < 15
            out = 55.03*exp(-0.13*15/(1+(15/273)))
        end
        return(out)
    end
    # Reaction norms
    WL     = CSV.read(raw"WL_re.csv", header=true, DataFrame)
    A_Long = CSV.read(raw"AMgam.csv", header=true, DataFrame)
    L_Surv = CSV.read(raw"Fin_LSurv.csv", header=true, DataFrame)
    L_Dur  = CSV.read(raw"LDgam.csv", header=true, DataFrame)
    #Misc parameters
    nTemp = 200; nDens = 200; nWing = 64;    # Number of environmental classes
    n     = nTemp * nDens;                   # Number of combinations of temp and density
    vert = 0.00                              # Latitude of simulated location
    ooo  = 16                                # Time taken to evaporate
    H_width  = 114                           # Width of container
    H_length = 114                           # Length of container
    height   = 38                            # Height of container
    surface_area = H_width * H_length        # Surface area of container
    m_precip     = height * surface_area     # Volume of container
    # Function that calculates moving averages
    function movingaverage(X::Vector,numofele::Int)
        BackDelta = div(numofele,2)
        ForwardDelta = isodd(numofele) ? div(numofele,2) : div(numofele,2) - 1
        len = length(X)
        Y = similar(X)
        for n = 1:len
            lo = max(1,n - BackDelta)
            hi = min(len,n + ForwardDelta)
            Y[n] = mean(X[lo:hi])
        end
        return Y
    end
    # Egg Q-class
    function Qui(t)
        if Water(t) < Water(t-1)
            out = 1 .* (m_precip   - Water(t))./(m_precip)
            if out > 1
                out = 1
            end
            if out < 0
                out = 0
            end
        else
            out = 0
        end
        return(out)
    end
    # Function for the temperature, smoothed
    function Temp(t )
        time = t   ;
        out = Temps_spl(time)
        return(out)
    end
    # Palette selection
    pal = get_color_palette(:auto, plot_color(:white));
    # Latitude: Guangzhou centre
    lat = latitude
    # Altitude (doesn't do anything anymore)
    alt  = 5
    # Number of days to be simulated
    tspan = (0.0, num_days);
    ttt = Int(tspan[2])                         #Number of days in simulation
    Water_Vect = repeat([0.0],ttt)              #Initialises precipitation vector
    modi = collect(range(0,1,length = ooo))     #Initialises vector for calculation
    Water_Vect[1:ooo] = repeat([m_precip],ooo);  #Starting conditions
    # Function that calculates precipitation accumulating in habitat at time t
    function Precip_check(t)
        out = Precip_spl.(t) .* surface_area
        if out < 0
            out = 0
        end
        return(out)
    end
    # Determines volume of water in habitat, accounting for evaporation and accumlation
    for i in (1 + ooo):Int(tspan[2])
        if Water_Vect[i-1] > 0
            Water_Vect[i] = Water_Vect[i-1] +   Precip_check(i) .+    EV_spl(i) .* surface_area
        else
            Water_Vect[i] = Water_Vect[i-1] + Precip_check(i)
        end
        if Water_Vect[i] .> m_precip
            Water_Vect[i] = m_precip
        end
    end
    Water_Vect[1:ooo] = repeat([m_precip],ooo);
    #Defines spline that calculates volume of water in habitat
    Water_spl = Spline1D(collect(range(1,Int(tspan[2]),length = Int(tspan[2]))),Water_Vect)
    # Function that calculates volume of water in habitat
    function Water(t)
        out =Water_spl.(  t)
        if out <0.1
            out = 0.1
        end
        return(( out ./1 ))
    end
    function resp(TEMP,time)
        10 .^ (0.45 + 0.095 * TEMP)   .* Water.(time) .* 1e-6 .+ 80
    end
    # Function that releases eggs from quiescence upon inundation
    function R(t)
        if Water(t) > Water(t-1)
            out  = 1 -  (m_precip ./1  - Water(t))./(m_precip ./ 1)
        else
            out = 0
        end
        if out < 0
            out = 0
        end
        if out > 1
            out = 1
        end
        return(out)
    end
    # Intiates vector for storing derivatives and adult recruitment
    derivs   = zeros(3 * nWing + 14)
    R_A      = zeros(nWing)
    R_A_I    = zeros(nWing);
    # Defines numbers for referencing derivatives
    numb_E_1   = Int(1)
    numb_E_D   = 2
    numb_E_Q   = 3
    numb_L_1   = 4
    numb_t_E   = 5
    numb_P_E   = 6
    numb_t_L   = 7
    numb_P_L   = 8
    numb_t_P   = 9
    numb_P_P   = 10
    numb_A_1   = 11
    numb_A_2       = nWing + 10
    numb_I_1       = nWing + 11
    numb_I_2       = 2* nWing + 10
    numb_EIP       = 2 * nWing + 11
    numb_P_EIP_1   = 2 * nWing + 12
    numb_P_EIP_2   = 3 * nWing + 11
    numb_H_S       = 3 * nWing + 12
    numb_H_I       = 3 * nWing + 13
    numb_H_R       = 3 * nWing + 14;
    # Creat300iscretisastion of environmental cues
    Temp_Disc          = collect(range(-15,40,length = nTemp+1))
    Temp_Vals          = (Temp_Disc[1:nTemp].+Temp_Disc[2:nTemp+1])./2
    Temp_Disc[1]       = -100
    Temp_Disc[nTemp+1] =  100
    Dens_Disc          = (collect(range(-8,1,length = nDens+1)))
    Dens_Vals          = ((Dens_Disc[1:nDens]+Dens_Disc[2:nDens+1])/2)
    Dens_Disc[nDens+1] = exp(14)
    Dens_Vals[1]       = -8
    Dens_Disc[1]       = -exp(14)
    Wing_Disc          = collect(range(1.5,4,length = nWing+1))
    Wing_Vals          = (Wing_Disc[1:nWing].+Wing_Disc[2:nWing+1])./2
    Wing_Disc[1]       = 0
    Wing_Disc[nWing+1] =  4;

    # Adult parameters
    # Defines adult wing lengths based on larval experience of temperature and density
    T_Vals = unique(WL[:,3])
    D_Vals = unique(WL[:,2])
    Wing_W       = transpose(reshape(WL[:,4], (400,400))); nodes = (T_Vals,D_Vals); sGird = interpolate(nodes,Wing_W,Gridded(Linear()));
    wing_spl     = extrapolate(sGird,Flat());
    Wing_Lengths = wing_spl(Temp_Vals,Dens_Vals); Wing_Lengths = reshape(transpose(Wing_Lengths), (1,n))[1,:];
    function wing_func(Temp,Dens)
        ifelse( Temp > 15, wing_spl(Temp,Dens),  wing_spl(15,Dens))
    end
    # Defines the mortality of adults based on wing length and adult temperature
    T_Vals = unique(A_Long[:,2])
    W_Vals = unique(A_Long[:,3])
    Long    = transpose(reshape(A_Long[:,4], (400,400))); nodes = (W_Vals,T_Vals); sGird = interpolate(nodes,Long,Gridded(Linear()));
    del_spl = extrapolate(sGird,Flat());
    function delta_AT(Temp)
        del = del_spl.(Wing_Vals,Temp)
        del =ifelse.(del .< 0.00001,0.00001,del)
        out = -log(0.5) ./ del
        out = ifelse.(out .> 0.5,  0.5, out)
        return(out)
    end
    #Defines the length of the gonotrophic cycle
    function G_Time(Temp)
        if Temp > 38.3
            Temp = 38.3
        end
        if Temp < 12
            Temp =  12
        end
        out = 1 ./ (1.93e-04  * Temp * (Temp - 10.25 ) * (38.32 - Temp) ^ (1/2))
    end
    #Defines function for eggs per gonotrophic cycle based on wing length
    function q_func(Temp)
        out =  exp.(2.35 .+ 0.69 .* Wing_Vals)./ 2
        return(out)
    end

    # Egg parameters
    # Egg development rate
    function g_E(Temp)
        if Temp < 6
            out =    -0.0008256     * 6 ^2  +0.0334072    * 6-0.0557825
        elseif Temp > 38
            out =    -0.0008256     * 38 ^2  +0.0334072    * 38 -0.0557825
        else
            out =    -0.0008256     * Temp ^2  +0.0334072    *Temp -0.0557825
        end
        return(out)
    end
    # Length of egg stage
    function tau_E(Temp)
        out = 1/g_E(Temp)
        return(out)
    end
    # Through egg stage survival
    function P_E(Temp)
        out = 12.217 * (1/(6.115*(2*pi)^0.5)) * exp(-0.5 * ((Temp - 24.672 )/6.115) ^ 2)
        if Temp < 6
            out = 12.217 * (1/(6.115*(2*pi)^0.5)) * exp(-0.5 * ((6 - 24.672 )/6.115) ^ 2)
        end
        if Temp > 40
            out = 12.217 * (1/(6.115*(2*pi)^0.5)) * exp(-0.5 * ((40 - 24.672 )/6.115) ^ 2)
        end
        return(out)
    end
    # Active egg mortality rate
    function delta_E(Temp)
        del = -log.(P_E(Temp)) ./ tau_E(Temp)
        if del > 0.99
            del = 0.99
        end
        return(del)
    end
    # Quiescent egg mortality rate
    function delta_E_Q(Temp)
        out = 12.217 * (1/(6.115*(2*pi)^0.5)) * exp(-0.5 * ((Temp - 24.672 )/6.115) ^ 2)
        if out <= 0.01
            out = 0.01
        end
        del = -log.(out) ./ tau_E(Temp)
        if del > 1
            del = 1
        end
        return(del)
    end
    # Diapausing egg mortality rate
    function delta_E_D(Temp)
        if Temp > -12
            out = 0.01
        else
            out = 0.1
        end
        return(out)
    end

    # Pupal parameters
    # Pupal development rate
    function g_P(TEMP)
        if TEMP > 40
            out = 2.916e-05  * 40 * (40- 1.008e+01 ) * (4.768e+01 - 40) ^ (1/8.317e-01 )
        elseif TEMP < 12.3
            out = 2.916e-05  *12.3 * (12.3- 1.008e+01 ) * (4.768e+01- 12.3) ^ (1/8.317e-01 )
        else
            out = 2.916e-05  * TEMP * (TEMP - 1.008e+01 ) * (4.768e+01 - TEMP) ^ (1/8.317e-01 )
        end
        return(out)
    end
    # Pupal development time
    function tau_P(TEMP)
        tau = 1/g_P(TEMP)
        return(tau)
    end
    # Through pupal-stage survival
    function P_P(TEMP)
        out =     -0.0070628      * TEMP ^2  +0.3331028    *TEMP -2.9878761
        if TEMP > 34
            out =  -0.0070628      *  34 ^2  +0.3331028    * 34 -2.9878761
        elseif TEMP < 12.3
            out =  -0.0070628      *  12.3 ^2  +0.3331028    * 12.3 -2.9878761
        end
        if out <= 0.000001
            out = 0.000001
        end
        if out > 0.99
            out = 0.99
        end
        return(out)
    end
    # Pupal mortality rate
    function delta_P(TEMP,time)
        del =   (- log.(P_P.(TEMP))./tau_P.(TEMP))
        if Precip_check.(time) .>=  m_precip  && del < 0.2 && Water(time) == m_precip
            del = 0.2
        end
        if del > 0.99
            del = 0.99
        end
        if Water(time) <= 6000
            del = 0.5
        end
        return(del)
    end

    # Larval parameters
    # Larval development time
    D_Vals = unique(L_Dur[:,2])
    T_Vals = unique(L_Dur[:,3])
    L_D       = transpose(reshape(L_Dur[:,4], (400,400))); nodes = (T_Vals,D_Vals); sGird = interpolate(nodes,L_D, Gridded(Linear()));
    tau_L_spl = extrapolate(sGird,Flat())
    # Larval development rate
    function g_L(TEMP,DENS)
        out = 1 ./tau_L_spl(TEMP,DENS)
        if TEMP < 12
            out = 1 ./ tau_L_spl(12, DENS)
        end
        return(out)
    end
    # Through larval-stage survival
    D_Vals = unique(L_Surv[:,3])
    T_Vals = unique(L_Surv[:,2])
    L_S       = (reshape(L_Surv[:,4], (400,400)));
    nodes = (T_Vals,D_Vals); sGird = interpolate(nodes,L_S, Gridded(Linear()));
    surv_L = extrapolate(sGird,Flat())
    function P_L(Temp,Dens)
        out = surv_L(Temp,Dens)
        if out < 0.01
            out = 0.01
        end
        if out > 0.99
            out = 0.99
        end
        return(out)
    end
    # Larval mortality rate
    # The biotic component of larval mortality, from temperature and copmetition
    function delta_L_spl(Temp,Dens)
        s = P_L(Temp,Dens)
        del = (- log.(s)./tau_L_spl(Temp,Dens))
        if del > 0.99
            del = 0.99
        end
        return(del)
    end
    # Modification of the larval mortality rate to account for abiotic factors such as overspill, evaporation, and overcrowding
    function delta_L(dens,temp,larv,time)
        mort = delta_L_spl(temp,dens)
        if Water(time) <= 6000 #evaporation
            mort = 0.5
        end
        out = mort .+    exp.( -   exp.( 1.0 .-  0.33.* abs.(larv + 1) ./ (Water.(time) ./ 1000))) #overcrowding
        if Precip_check.(time) .>=   m_precip  && out < 0.2 && Water(time) == m_precip #overspill induced flushing
            out = 0.2
        end
        if out > 0.99
            out = 0.99
        end
        if out < 0.01
            out = 0.01
        end
        return(out)
    end
    # Inoculation function
    function on_off(t,a)
        impulse = zeros(length(a));
        for i = 1:length(a)
            if t <= (a[i] .+ 10) && t .>= (a[i] .- 10)
                impulse[i] = 1
            end
        end
        return(impulse)
    end
    # Function for assigning adults to environmental classes
    function w_func(TEMP,DENS,VECT)
        VECT_TEMP = VECT
        inp  = wing_func(TEMP,DENS)
        for z in 1:nWing
            if Wing_Disc[z] <= inp < Wing_Disc[z + 1]
                VECT_TEMP[z] = 1
            end
        end
        return(VECT_TEMP)
    end

    # Function that calculates photoperiod at time t
    function PP(T)
        L = lat
        EPS = asin(0.39795 * cos(0.2163108 + 2 * atan(0.9671396 * tan(0.00860 * (T-183.5)))))
        NUM = sin(0.8333 *pi/180) + (sin(L*pi/180) * sin(EPS))
        DEN = cos(L*pi/180) * cos(EPS)
        DAYLIGHT = 24 - (24/pi) * acos(NUM / DEN)
        return(DAYLIGHT)
    end
    # Defines critical photoperiod
    a = 1.5 * abs(lat)/14 + 12.5 - 28*1.5/14;
    # Function for production of diapausing eggs
    function Dia(t)
        if PP.(t) < PP(t-1)
            if Temp(t) < 18
                out =  1- 1/(1 + 15*exp( (PP(t) - a)))
            else
                out = 1
            end
        else
            out = 1
        end
        return(out)
    end
    # Function for release from diapause
    function H(t)
        if Temp(t) > 14.5
            if PP(t) > a
                if PP(t) > PP(t-1)
                    out = 1
                else
                    out = 0
                end
            else
                out = 0
            end
        else
            out = 0
        end
    end
    
    # Seeding function
    seed_num = Dates.value(Date(seed_date, "dd.mm.yyyy") - Date(t_start, "dd.mm.yyyy"))
    ## Function for initial impulse
    function infect(time)
        # Convert simulation time into the day of the year
        #day_of_year = mod(time, 365)
        impulse = 0  # Default value

        # seeding 
    
        if time >= seed_num - 1 && time <= seed_num 
            impulse = 1
        end
    
        return impulse
    end

    #Defines intial conditions
    # Fano susceptible population size (tested 12000)
    pop = 7000;
    y0                             = zeros(3 * nWing + 14)
    y0[numb_E_1]                   = 1
    y0[numb_E_D]                   = 150
    y0[numb_E_Q]                   = 1
    y0[numb_L_1]                   = 100
    y0[numb_t_P]                   = tau_P(Temp(0))
    y0[numb_t_L]                   = tau_L_spl(Temp(0),log(Water(0)./y0[numb_L_1]))
    y0[numb_t_E]                   = tau_E(Temp(0))
    y0[numb_P_P]                   = exp(-y0[numb_t_P] * delta_P(Temp(0),0))
    y0[numb_P_L]                   = exp(-y0[numb_t_L] * delta_L(log(Water(0)./y0[numb_L_1]),Temp(0),y0[numb_L_1],0))
    y0[numb_P_E]                   = exp(-y0[numb_t_E] * delta_E(Temp(0)))
    y0[numb_EIP]                   = EIP(Temp(0))
    y0[numb_P_EIP_1:numb_P_EIP_2] .= exp.(-y0[numb_EIP] .* delta_AT(Temp(0)))
    y0[numb_H_S]                   = pop
    y0[numb_H_I]                   = 0
    lags                           = y0
    w_A                            = zeros(nWing)
    h(p, t)                        = lags
    R_A_store                      = [];

    # incubation (iip) and recovery time
    t_iip = 3
    t_recov = 7
    
    # Function that defines the model
    function aedes_model(dy,y,h,p,t)
        #Temperatures at developmental milestones
        TEMP_NOW = Temp(t              )[1] ;
        TEMP_E   = Temp(t - y[numb_t_E])[1] ;
        TEMP_L   = Temp(t - y[numb_t_L])[1] ;
        TEMP_P   = Temp(t - y[numb_t_P])[1] ;
        TEMP_EIP = Temp(t - y[numb_EIP])[1] ;
        #State of variables at developmental milestones
        ylag_E   = h(p,t -  y[numb_t_E])    ;
        ylag_L   = h(p,t -  y[numb_t_L])    ;
        ylag_P   = h(p,t -  y[numb_t_P])    ;
        ylag_EIP = h(p,t -  y[numb_EIP])    ;
        ylag_inf = h(p,t - t_iip)    ;
        ylag_rev = h(p,t - t_iip - t_recov);
        #Average temperature over the larval period
        t1 = t-y[numb_t_P]
        t2 = t-y[numb_t_P]-ylag_P[numb_t_L]
        TEMP_P_AVG = mean(Temp.(range(t2,t1 , length = 20)))
        #Temperatures at developmental milestones
        TEMP_E_L   = Temp(t - y[numb_t_L] - ylag_L[numb_t_E])[1] ;
        TEMP_L_P   = Temp(t - y[numb_t_P] - ylag_P[numb_t_L])[1] ;
        #State of variables at developmental milestones
        ylag_E_L   = h(p, t - y[numb_t_L] - ylag_L[numb_t_E])    ;
        ylag_L_P   = h(p, t - y[numb_t_P] - ylag_P[numb_t_L])    ;
        ylag_L_L   = h(p, t - y[numb_t_L] - ylag_L[numb_t_L])    ;
        #Temperatures at developmental milestones
        TEMP_E_L_P = Temp(t - y[numb_t_P] - ylag_P[numb_t_L] - ylag_L_P[numb_t_E] )[1] ;
        #State of variables at developmental milestones
        ylag_E_L_P =  h(p,t - y[numb_t_P] - ylag_P[numb_t_L] - ylag_L_P[numb_t_E] )    ;
        #Calculates integrals for the intensity of larval competition over the larval period at various times
        v         = log(             ( resp(TEMP_NOW,t)) / abs(       y[numb_L_1] + 1)) ;
        v_L       = log(             (  resp(TEMP_L,t - y[numb_t_L])  ) / abs(  ylag_L[numb_L_1] + 1)) ;
        v_L_P     = log(             (  resp(TEMP_L_P, t - ylag_P[numb_t_L] - y[numb_t_P])) / abs(ylag_L_P[numb_L_1] + 1)) ;
        v_P       = log(             (  resp(TEMP_P, t - y[numb_t_P])  ) / abs(  ylag_P[numb_L_1] + 1)) ;
        val_P     = log( quadgk(t -> (  resp(Temp(t),t) ) / abs(  h(p,t)[numb_L_1] + 1)  , t - y[numb_t_P] - ylag_P[numb_t_L], t - y[numb_t_P], rtol = 1)[1] ./ ylag_P[numb_t_L] ) ;
        #Proportion of adults producing diapausing eggs
        Dia_Now   = Dia(t                                                      ) ;
        Dia_E     = Dia(t - y[numb_t_E]                                        ) ;
        Dia_E_L   = Dia(t - y[numb_t_L] - ylag_L[numb_t_E]                     ) ;
        Dia_E_L_P = Dia(t - y[numb_t_P] - ylag_P[numb_t_L] - ylag_L_P[numb_t_E]) ;
        #Defines transition functions
        w_A_temp  = zeros(nWing)                              ;
        w_A       = w_func(TEMP_P_AVG, val_P, w_A_temp) ;
        #Recruitement terms
        R_E       = Dia_Now .* sum(repeat(q_func(TEMP_NOW),2) .* y[numb_A_1:numb_I_2] ./ (G_Time(TEMP_NOW)))  .+ on_off(t,0)[1] ;
        R_E_lag   = Dia_E   .* ( g_E(TEMP_NOW) / g_E(TEMP_E) ) .* sum( repeat(q_func(TEMP_E),2) .* ylag_E[numb_A_1:numb_I_2] ./ (G_Time(TEMP_E)) ) .* y[numb_P_E] .+
        ( g_E(TEMP_NOW) / g_E(TEMP_E) ) .* on_off(t-y[numb_t_E],0)[1] .* y[numb_P_E] ;
        R_E_D     =   (1-Dia_Now) .* (sum(repeat(q_func(TEMP_NOW),2) .* y[numb_A_1:numb_I_2] ./ (G_Time(TEMP_NOW)))) ;
        R_E_D_lag =  H(t) * y[numb_E_D] ;
        R_E_Q     = Qui(t) .* (R_E_lag .+ R_E_D_lag)
        R_E_Q_mat =  R(t) * y[numb_E_Q]
        R_L     = (1-Qui(t )) * Dia_E .* (g_E(TEMP_NOW) /  g_E(TEMP_E) ) * sum(repeat(q_func(TEMP_E),2) .* ylag_E[numb_A_1:numb_I_2] ./ (G_Time(TEMP_E))) .* y[numb_P_E]  .+
        (1-Qui(t )) *  H(t) * y[numb_E_D] .+
        R(t) * y[numb_E_Q] .+
        (g_E(TEMP_NOW) /  g_E(TEMP_E) ) * on_off(t-y[numb_t_E],0)[1] .* y[numb_P_E]
        R_L_lag = (1-Qui(t  - y[numb_t_L])) * Dia_E_L .* (g_E(TEMP_L ) / g_E(TEMP_E_L)) * (g_L(TEMP_NOW, v ) ./ g_L(TEMP_L, v_L )) * sum(repeat(q_func(TEMP_E_L),2) .* ylag_E_L[numb_A_1:numb_I_2] ./ (G_Time(TEMP_E_L))) .* ylag_L[numb_P_E] .* y[numb_P_L] .+
        (1-Qui(t  - y[numb_t_L])) *(g_L(TEMP_NOW, v ) ./ g_L(TEMP_L, v_L )) * H(t - y[numb_t_L]) * ylag_L[numb_E_D] .* y[numb_P_L] .+
        (g_L(TEMP_NOW, v ) ./ g_L(TEMP_L, v_L )) * R(t - y[numb_t_L]) * ylag_L[numb_E_Q] .* y[numb_P_L] .+
        (g_E(TEMP_L ) / g_E(TEMP_E_L)) * (g_L(TEMP_NOW, v ) ./ g_L(TEMP_L, v_L )) * on_off(t-ylag_L[numb_t_E] - y[numb_t_L],0)[1].*  ylag_L[numb_P_E]  * y[numb_P_L]
        R_A   = (1-Qui(t  - ylag_P[numb_t_L] - y[numb_t_P])) * Dia_E_L_P .* w_A .* (g_E(TEMP_L_P)/ g_E(TEMP_E_L_P)) * (g_L(TEMP_P, v_P) ./ g_L(TEMP_L_P, v_L_P)) * (g_P(TEMP_NOW) ./ g_P(TEMP_P)) .*sum(repeat(q_func(TEMP_E_L_P),2) .*   ylag_E_L_P[numb_A_1:numb_I_2] ./ (G_Time(TEMP_E_L_P))) .* ylag_L_P[numb_P_E] .* ylag_P[numb_P_L] .* y[numb_P_P] .+
        (1-Qui(t  - ylag_P[numb_t_L] - y[numb_t_P]))*w_A .*(g_L(TEMP_P, v_P) ./ g_L(TEMP_L_P, v_L_P)) * (g_P(TEMP_NOW) ./ g_P(TEMP_P)) .*H(t - ylag_P[numb_t_L] - y[numb_t_P]) .* ylag_L_P[numb_E_D] .* ylag_P[numb_P_L] .* y[numb_P_P] .+
        w_A .*(g_L(TEMP_P, v_P) ./ g_L(TEMP_L_P, v_L_P)) * (g_P(TEMP_NOW) ./ g_P(TEMP_P)) .*R(t - ylag_P[numb_t_L] - y[numb_t_P]) .* ylag_L_P[numb_E_Q].* ylag_P[numb_P_L] .* y[numb_P_P] .+
        w_A .* (g_E(TEMP_L_P)/ g_E(TEMP_E_L_P)) * (g_L(TEMP_P, v_P) ./ g_L(TEMP_L_P, v_L_P)) * (g_P(TEMP_NOW) ./ g_P(TEMP_P)) .* on_off(t-ylag_L_P[numb_t_E] - ylag_P[numb_t_L] - y[numb_t_P],repeat([0],nWing)) .* ylag_L_P[numb_P_E] * ylag_P[numb_P_L] * y[numb_P_P]
        #Differential equations to be solved
        dy[numb_E_1]          = R_E     - R_E_lag     -  delta_E(TEMP_NOW) * y[numb_E_1]                 ;
        dy[numb_E_D]          = R_E_D   - R_E_D_lag   -  delta_E_D(TEMP_NOW) * y[numb_E_D]                 ;
        dy[numb_E_Q]          = R_E_Q   - R_E_Q_mat   -  delta_E_Q(TEMP_NOW) * y[numb_E_Q]                 ;
        dy[numb_L_1]          = R_L   - R_L_lag   -  delta_L(v,TEMP_NOW,y[numb_L_1],t)  * y[numb_L_1]          ;
        dy[numb_A_1:numb_A_2] = R_A              .-  delta_AT(TEMP_NOW )         .* y[numb_A_1:numb_A_2] .- ((1/EIP(TEMP_NOW)) ./ (1/EIP(TEMP_EIP))) * (human_vector(TEMP_EIP) .* bite(TEMP_EIP)) .* ylag_EIP[numb_A_1:numb_A_2]  .*  (ylag_EIP[numb_H_I]  ./ (ylag_EIP[numb_H_S]   + ylag_EIP[numb_H_I]  + ylag_EIP[numb_H_R] .+ 15000   )) .* y[numb_P_EIP_1:numb_P_EIP_2]                       ;
        dy[numb_I_1:numb_I_2] =                  .-  delta_AT(TEMP_NOW )         .* y[numb_I_1:numb_I_2] .+ ((1/EIP(TEMP_NOW)) ./ (1/EIP(TEMP_EIP))) * (human_vector(TEMP_EIP) .* bite(TEMP_EIP)) .* ylag_EIP[numb_A_1:numb_A_2]  .*  (ylag_EIP[numb_H_I]  ./ (ylag_EIP[numb_H_S]   + ylag_EIP[numb_H_I]  + ylag_EIP[numb_H_R] .+ 15000   )) .* y[numb_P_EIP_1:numb_P_EIP_2] ;
        dy[numb_t_E]          =  1 - (g_E(TEMP_NOW        ) ./ g_E(TEMP_E)     )      ;
        dy[numb_t_L]          =  1 - (g_L(TEMP_NOW,     v ) ./ g_L(TEMP_L, v_L))      ;
        dy[numb_t_P]          =  1 - (g_P(TEMP_NOW        ) ./ g_P(TEMP_P)     )      ;
        dy[numb_EIP]          =  1 - ((1/EIP(TEMP_NOW)) ./ (1/EIP(TEMP_EIP)))
        dy[numb_P_E]                  = y[numb_P_E]   * ((g_E(TEMP_NOW   ) * delta_E(TEMP_E))                       / g_E(TEMP_E)     - delta_E(TEMP_NOW)               ) ;
        dy[numb_P_L]                  = y[numb_P_L]   * ((g_L(TEMP_NOW, v) * delta_L(v_L,TEMP_L, ylag_L[numb_L_1], t - y[numb_t_L])) / g_L(TEMP_L,v_L) - delta_L(v,TEMP_NOW, y[numb_L_1],t)) ;
        dy[numb_P_P]                  = y[numb_P_P]   * ((g_P(TEMP_NOW   ) * delta_P(TEMP_P, t - y[numb_t_P]))                       / g_P(TEMP_P)     - delta_P(TEMP_NOW,t)               ) ;
        dy[numb_P_EIP_1:numb_P_EIP_2] = y[numb_P_EIP_1:numb_P_EIP_2] .* (((1/EIP(TEMP_NOW))    * delta_AT(TEMP_EIP))    ./ (1/EIP(TEMP_EIP))  .- delta_AT(TEMP_NOW)              ) ;
        dy[numb_H_S]                  = .- vector_human(TEMP_NOW) .* bite(TEMP_NOW)         .*  40* 4  .* abs(sum(       y[numb_I_1:numb_I_2])) .*        y[numb_H_S] ./ (y[numb_H_S]        + y[numb_H_I]        + y[numb_H_R]    .+ 15000       )
        dy[numb_H_I]                  = .+ vector_human(Temp(t - t_iip)) .* bite(Temp(t - t_iip))   .*  40* 4  .* abs(sum(ylag_inf[numb_I_1:numb_I_2])) .* ylag_inf[numb_H_S] ./ (ylag_inf[numb_H_S] + ylag_inf[numb_H_I] + ylag_inf[numb_H_R]  .+ 15000     ) .+ infect(t) .- infect(t-t_recov)  .-  (vector_human(Temp(t - t_iip - t_recov)) .* bite(Temp(t - t_iip - t_recov))) .*   40* 4  .* sum( ylag_rev[numb_I_1:numb_I_2]) .* ylag_rev[numb_H_S] ./ (ylag_rev[numb_H_S] + ylag_rev[numb_H_I] + ylag_rev[numb_H_R]  .+ 15000     )
        dy[numb_H_R]                  = .+ vector_human(Temp(t - t_iip - t_recov)) .* bite(Temp(t - t_iip - t_recov))   .*  40* 4  .* abs(sum(ylag_rev[numb_I_1:numb_I_2])) .* ylag_rev[numb_H_S] ./ (ylag_rev[numb_H_S] + ylag_rev[numb_H_I] + ylag_rev[numb_H_R]  .+ 15000      ) .+ infect(t-t_recov)
        #push!(R_A_store,R_A)
    end
    #Defines named parameters
    p = ( n,Water, tau_E, tau_P, delta_E, delta_L, P_E, numb_L_1, numb_A_1, numb_A_2, nDens)
    prob = DDEProblem(aedes_model,y0,h,tspan,p)
    # Remark: Sandeep-> lag declaration is removed as it is not required in RK4() Algorithm.
    alg = MethodOfSteps(RK4())
    sol = solve(prob, alg, reltol = 1e-3, abstol = 1e-3, saveat = 1.0)

    thing = sol[numb_A_1:numb_A_2,:]
    Temp_sim = Temp.(sol.t)
    Temp_sim = ifelse.(Temp_sim .> 37 ,37 ,Temp_sim)
    Delta_sim = zeros(nWing,length(sol.t))
    for i = 1:length(sol.t)
        Delta_sim[:,i] .= delta_AT(Temp_sim[i])
    end

    R_0 = reshape(repeat(Float64[0], Int64(tspan[2]) * nWing) , (nWing,Int64(tspan[2])))
    ADU = reshape(repeat(Float64[0], Int64(tspan[2]) * nWing) , (nWing,Int64(tspan[2])))
    TIME_R_0 = Spline1D(sol.t, sol.t .- sol[numb_EIP,:])
    TAU_EIP = Spline1D(sol.t, sol[numb_EIP,:])
    for j in 1:nWing
        store_ad = Spline1D(sol.t, thing[j,:])
        store_del = Spline1D(sol.t, 1 ./ Delta_sim[j,:])
        store_surv = Spline1D(sol.t, sol[numb_P_EIP_1 + j - 1,:])
        #store_H_S = Spline1D(sol.t, sol[numb_H_S,:])
        #store_H_I = Spline1D(sol.t, sol[numb_H_I,:])
        #store_H_R = Spline1D(sol.t, sol[numb_H_R,:])
        #store_A_S = Spline1D(sol.t, sum(thing .* sol[numb_P_EIP_2:numb_P_EIP_2,:], dims = 1)[1,:])
        times = collect(range(1,tspan[2],length = Int64(tspan[2])))
        for i in 1:Int64(tspan[2])
            low = Float64(times[i] )
            high = Float64(times[i] + 4)
            R_0[j,i] = quadgk(t -> 40 * 4 * ((1/EIP(Temp(t))) ./ (1/EIP(Temp(t - TAU_EIP(t))))) * human_vector(Temp(t - TAU_EIP(t))) * bite(Temp(t - TAU_EIP(t))) * store_ad(t - TAU_EIP(t)) * store_surv(t) * (1 ./  (y0[numb_H_S]  + 15000 ))
            *  quadgk(u -> bite(Temp(u)) * vector_human(Temp(u)) * (y0[numb_H_S] ./  (y0[numb_H_S]  + 15000 )),  t, t +  store_del(t), rtol = 1)[1], low, high, rtol = 1)[1]
            ADU[j,i] = store_ad(i)
        end
    end
    R_0_spl = Spline1D(TIME_R_0.(collect(range(1,Int(tspan[2]),length = Int(tspan[2])))), sum(R_0',dims = 2)[:,1]);

    # sum of adult vectors
    sumA = sum(sol[j,:] for j = numb_A_1:numb_A_2)

    # Extrinsic incubation period
    tauE = sol[numb_t_E, :]

    # Computing daily new number of cases (number of cases on the day-t: I[t]-I[t-1])
    daily_new_cases = [0.0]
    for t in 2:length(sol.t)
        new_cases = sol[numb_H_I,:][t]-sol[numb_H_I,:][t-1]
        push!(daily_new_cases, 5*max(new_cases, 0.0))
    end
    
    # writing output time series
    start_date = Date(t_start, "dd.mm.yyyy")
    dates =  start_date .+ Day.(Int.(sol.t))
    infectious = sol[numb_H_I,:]
    recovered = sol[numb_H_R,:]
    rt =  R_0_spl.(sol.t)
    sim_data = hcat(sumA, infectious, daily_new_cases, recovered, tauE, rt)
    ts = TimeArray(dates, sim_data, ["m_total", "i_inst", "i_daily", "Hr", "tE", "Rt"])
    return ts 

end