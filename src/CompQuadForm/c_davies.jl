function log1(x::Float64, first::Bool) 
    if abs(x) > 0.1
        if first == true
            return(log(1.0 + x))
        else 
            return(log(1.0 +x) - x)
        end
    else 
        y = x / (2.0 + x)
        term = 2.0 * y^3
        k = 3.0
        s = (first ? 2.0 : - x) * y
        y = y^2
        for s1 = s + term / k
            if s1 != s
                k = k + 2.0
                term = term * y
                s = s1
            end
        end
        return(s)
    end
end

#To avoid underflows
function exp1(x::Float64)               
    if x < -50.0
        return(0)
    else 
        return(exp(x))
    end
end

#find bound on tail probability using mgf, cutoff point returned to cx
function errbd!(u::Float64, cx::Vector, sigsq::Float64, n::Vector, lb::Vector, nc::Vector, r::Int64) 
    if size(n, 1) != r
        error("size of n does not equal r")
        end
    xconst = u * sigsq 
    sum1 = u * xconst
    u = 2.0 * u
    for j in r:1
        nj = n[j]     
        lj = lb[j] 
        ncj = nc[j] 
        x = u * lj 
        y = 1.0 - x  
        xconst = xconst + lj * (ncj / y + nj) / y
        sum1 = sum1 + ncj * (x/y)^2 + nj * x^2 /y +log1(-x, false) 
    end
    xconst = cx[1]
    return(exp1(-0.5 * sum1))
end

#find ctff so that p(qf > ctff) < accx  if (upn > 0, p(qf < ctff) < accx otherwise
function ctff(accx::Float64, upn::Float64, meanvalue::Float64, lmax::Float64, lmin::Float64, sigsq::Float64, n::Vector, lb::Vector, nc::Vector, r::Int64, xconst::Int64) 
    u1 = 0.0
    c1 = meanvalue 
    c2 = Number[]
    push!(c2, 0)
    xconst = [xconst]
    if upn > 0.0
        rb = 2.0 * lmax 
    else
        rb = 2.0 * lmin 
    end
    for u = upn / (1.0 + upn * rb)
        while errbd!(u, c2, sigsq, n, lb, nc, r) > accx 
            u1 = upn
            c2[1] = xconst[1]
            c1 = c2[1]
            upn = 2.0 * upn
        end
    end
    for u = (c1 - meanvalue) / (c2[1] - meanvalue)
        while u < 0.9
            u = (u1 + upn) / 2.0
            if errbd!(u / (1.0 + u * rb), xconst, sigsq, n, lb, nc, r) > accx 
                u1 = u
                c1 = xconst[1]
             else
                upn = u
                c2[1] = xconst[1]  
             end
        end 
    end
    ctff = c2[1] 
    return(ctff)
end

#bound integration error due to truncation at u
function truncation(u::Float64, tausq::Float64, lb::Vector, nc::Vector, n::Vector, sigsq::Float64, r::Int64) 
    if size(lb, 1) != r
        error("size of lb does not equal r")
    elseif size(nc, 1) != r
        error("size of nc does not equal r")
    elseif size(n, 1) != r
        error("size of n does not equal r")
    end
    sum1 = 0.0
    prod2 = 0.0
    prod3 = 0.0
    s = 0.0
    sum2 = (sigsq +tausq) * u^2
    prod1 = 2.0 * sum2
    u = 2.0 * u
    for j in 1:r 
        lj = lb[j] 
        ncj = nc[j]
        nj = n[j]
        x = (u * lj)^2
        sum1 = sum1 + ncj * x / (1.0 + x)
        if x > 1.0
            prod2 = prod2 + nj * log(x)
            prod3 = prod3 + nj * log1(x, true) 
            s = s + nj
        else
            prod1 = prod1 + nj * log1(x, true) 
        end
    end
    sum1 = 0.5 *sum1 
    prod2 = prod1 + prod2
    prod3 = prod1 + prod3
    x = exp1(-sum1 - 0.25 * prod2) / pi
    y = exp1(-sum1 - 0.25 * prod3) / pi
    if s == 0
        err1 = 1.0
    else
        err1 = x * 2.0 / s
    end
    err2 = (prod3 > 1.0) ? 2.5 * y : 1.0
    if err2 < err1 
        err1 = err2
        x = 0.5 * sum2
    end
    err2 =  ( x  <=  y )  ? 1.0  : y / x
    return(err1 < err2 ? err1 : err2)
end
#truncation(10.0/4.0, 0.0, [1, 2, 3], [1, 2, 3], [1, 2, 3], 7.0, 3)

#find u such that truncation(u) < accx and truncation(u / 1.2) > accx
function findu(ut::Float64, accx::Float64, tausq::Float64, lb::Vector, nc::Vector, n::Vector, sigsq::Float64, r::Int64) 
    divis = [2.0, 1.4, 1.2, 1.1]
    u = ut / 4
    if truncation(u, 0.0, lb, nc, n, sigsq, r) > accx
        for u = ut
            if truncation(u, 0.0, lb, nc, n, sigsq, r) > accx 
                ut = ut * 4.0
            end 
        end
    else
        ut = u
        for u = u / 4.0
            if truncation(u, 0.0, lb, nc, n, sigsq, r) <= accx 
                ut = u 
            end
        end
    end
    for i in 1:4
        u = ut / divis[i]
        if truncation(u, 0.0, lb, nc, n, sigsq, r) <= accx 
            ut = u
        end
    end
end
#findu(10.0, 9.0, 8.0, [1, 2, 3], [1, 2, 3], [1, 2, 3], 7.0, 3)

#carry out integration with nterm terms, at stepsize interv.  
#if (! mainx) multiply integrand by 1.0 - exp(-0.5 * tausq * u ^ 2)
function integrate(nterm::Int64, interv::Float64, tausq::Float64, mainx::Bool, c::Float64, sigsq::Float64, n::Vector, lb::Vector, nc::Vector, r::Int64) 
    if size(lb, 1) != r
        error("size of lb does not equal r")
    elseif size(nc, 1) != r
        error("size of nc does not equal r")
    elseif size(n, 1) != r
        error("size of n does not equal r")
    end
    inpi = interv / pi
    for k in nterm:0
        u = (k + 5.0) * interv
        sum1 = -2.0 * u * c 
        sum2 = abs(sum1)
        sum3 = - 0.5 * sigsq * u^2
        for j in r:1 
            nj = n[j]
            x = 2.0 * lb[j] * u
            y = x^2
            sum3 = sum3 - 0.25 * nj * log1(y, true) 
            y = nc[j] * x / (1.0 + y)
            z = nj * atan(x) + y
            sum1 = sum1 + z
            sum2 = sum2 + abs(z)
        end
        x = inpi * exp1(sum3) / u
        if mainx == false 
            x = x * (1.0 - exp1(-0.5 * tausq * u^2))
        end
        sum1 = sin(0.5 * sum1) * x
        sum2 = 0.5 * sum2 * x
        intl = 0.0 + sum1 
        ersm = 0.0 + sum2 
    end
end
#integrate(10, 9.0, 8.0, false, 7.0, 6.0, [1, 2, 3], [1, 2, 3], [1, 2, 3], 3)

#find order of absolute values of lb
function order(lb::Vector, th::Vector, r::Int64, j::Int64)
    if size(lb, 1) != r
        error("size of lb does not equal r")
    end
    for j in 1:r 
        lj = abs(lb[j])
        for k in j:1  
            if lj > abs(lb[th[k]]) 
                th[k + 1] = th[k]
            else
                th[k + 1] = j
            end
        k = -1
        end
    end
    ndstart = false
end
#order([1, 2, 3], [1, 2, 3], 3, 3)

#coef of tausq in error when convergence factor of
#exp1(-0.5 * tausq * u ^ 2) is used when df is evaluated at x
function cfe(x::Float64, lb::Vector, th::Vector, r::Int64, ndstart::Bool, n::Vector, nc::Vector, j::Int64) 
    #may need error for not equal here
    if ndstart == true
        order(lb, th, r, j) 
    end
    axl = abs(x) 
    sxl = (x > 0.0) ? 1.0 : -1.0 
    sum1 = 0
    for m in r:1
        t = th[m]
        if lb[t] * sxl > 0.0
            lj = abs(lb[t])
            axl1 = axl - lj * (n[t] + nc[t])
            axl2 = lj/(log(2.0)/8)
            if axl > axl2
                axl = axl2
                sum1 = (axl - axl1) / lj
            end
            for k in m:1
                sum1 = sum1 + (n[th[k]] + n[th[k]])
                @goto l
            end
        end
    end
    @label l
        if sum1 > 100.0 
            fail = true
            return(1.0)
        else
            return(2.0^((sum1 / 4.0) / (pi * axl^2)))
        end
end
#cfe(10.0, [1, 2, 3], [1, 2, 3], 3, false, [1, 2, 3], [1, 2, 3], 3)

#distribution function of a linear combination of non-central chi-squared random variables :
 #input:
    #lb[j]            coefficient of j-th chi-squared variable
    #nc[j]            non-centrality parameter
    #n[j]             degrees of freedom
    #j = 0, 2 ... r-1
    #sigma            coefficient of standard normal variable
    #c                point at which df is to be evaluated
    #lim              maximum number of terms in integration
    #acc              maximum error
 
 #output:
    #ifault = 1       required accuracy NOT achieved
             #2       round-off error possibly significant
             #3       invalid parameters
             #4       unable to locate integration parameters
             #5       out of memory
 
    #trace[0]         absolute sum
    #trace[1]         total number of integration terms
    #trace[2]         number of integrations
    #trace[3]         integration interval in final integration
    #trace[4]         truncation point in initial integration
    #trace[5]         s.d. of initial convergence factor
    #trace[6]         cycles to locate integration parameters
function qfc(lb::Vector, nc::Vector, n::Vector, s::Int64, r::Int64, sigma::Float64, c::Float64, lim::Int64, acc::Float64, ifault::Int64, res::Float64, tausq::Float64, th::Vector, ndstart::Bool) 
    lim = lim[1]
    c = c[1]
    xconst = Number[]
    xconst = 0.0

    #Label L1, L2
    trace = zeros(Float64, 7)
    ifault = 0 
    count = 0
    intl = 0.0
    ersm = 0.0
    qfval = -1.0
    acc1 = acc[1]
    ndtsrt = true
    fail = false
    xlim = lim
    
    # find mean, sd, max and min of lb,
    #check that parameter values are valid 
    sigsq = (sigma[1])^2 
    sd = sigsq
    lmax = 0.0
    lmin = 0.0
    meanvalue = 0.0
    for j in 1:r 
        nj = n[j]  
        lj = lb[j]
        ncj = nc[j]
        if nj < 0 || ncj < 0.0 
            ifault = 3
            @goto endofproc
        end
        sd = sd + (lj)^2 * ((2 * nj) + 4.0 * ncj)
        meanvalue = meanvalue + lj * (nj + ncj)
        if lmax < lj
            lmax = lj
        elseif lmin > lj 
            lmin = lj
        end
    end
    if sd == 0.0
        qval = c > 0.0 ? 1.0 : 0.0
        @goto endofproc
    end
    if lmin == 0.0 && lmax == 0.0 && sigma[1] == 0.0
        ifault = 3
        @goto endofproc
    end
    sd = sqrt(sd)
    almx = (lmax < -lmin) ? - lmin : lmax
    
    #starting values for findu, ctff
    utx = 16.0 / sd
    up = 4.5 / sd
    un = -up
    
    #truncation point with no convergence factor
    findu(utx, 0.5 * acc1, tausq, lb, nc, n, sigsq, r)
    
    #Does convergence factor help?
    if c != 0.0 && almx > 0.07 * sd
        tausq = 0.25 * acc1 / cfe(c, lb, th, r, ndstart, n, nc, s)
        if fail == true
            fail = false
        elseif truncation(utx, tausq, lb, nc, n, sigsq, r) < 0.2 * acc1
            sigsq = sigsq + tausq
            findu(utx, 0.25 * acc1, tausq, lb, nc, n, sigsq, r)
            trace[6] = sqrt(tausq) 
        end
    end
    trace[5] = utx 
    acc1 = 0.5 * acc1

    #find RANGE of distribution, quit if outside this
    @label l1
        d1 = ctff(acc1, up[1], meanvalue, convert(Float64, lmax), lmin, sigsq, n, lb, nc, r, convert(Int64, xconst)) - c    
        qfval = -1 
        if d1 < 0.0
            qfval = 1.0
            @goto endofproc
        end
        d2 = c - ctff(acc1, un[1], meanvalue, convert(Float64, lmax), lmin, sigsq, n, lb, nc, r, convert(Int64, xconst)) 
        if d2 < 0.0
            qfval = 0.0
            @goto endofproc
        end
        ##randomly stops in the above section may be because there is not a call to this label above


        #Find integration Interval
        intv = 2.0 * pi / ((d1 > d2) ? d1 : d2)
        
        #calculate number of terms required for main and auxillary integrations
        xnt = utx / intv
        xntm = 3.0 / sqrt(acc1)
        if xnt > xntm * 1.5
            
            #Parameters for auxillary integrations
            if xntm > xlim 
                ifault = 1 
                @goto endofproc
            end
            ntm = floor(xntm + 0.5) 
            intv1 = utx / ntm
            x = 2.0 * pi / intv1 
            if x <= abs(c) 
                @goto l2
            end
            
            #Calculate Convergence Factor
            tausq = 0.33 * acc1 / (1.1 * (cfe(c-x, lb, th, r, ndstart, n, nc) + cfe(c+x, lb, th, r, ndstart, n, nc)))
            if fail
                @goto l2

            end
            acc1 = 0.67 * acc1
    
            
            #auxillary integration
            integrate(convert(Int64, ntm), intv1, tausq, false, c, sigsq, n, lb, nc, r)
            xlim = xlim - xntm
            sigsq = sigsq + tausq
            trace[3] = trace[3] + 1 
            trace[2] = trace[2] + ntm + 1
            
            #find truncation point with new convergence factor
            findu(utx, 0.25 * acc1, tausq, lb, nc, n, sigsq, r)
            acc1 = 0.75 * acc1
            @goto l1
        end
    
        #Main Integration
    @label l2
        trace[4] = intv 
        if xnt > xlim
            ifault = 1
            @goto endofproc
        end
        nt = floor(xnt + 0.5) 
        integrate(convert(Int64, nt), intv, 0.0, true, c, sigsq, n, lb, nc, r)
        trace[3] = trace[3] + 1 
        trace[2] = trace[2] + nt + 1
        qfval = 0.5 - intl 
        trace[1] = ersm 
        
        #test whether round-off error could be significant 
        #allow for radix 8 or 16 machines
        up = ersm
        x = up + acc[1] / 10.0 
        rats = [1, 2, 4, 8]
        for j in 0:4
            if rats[j] * x == rats[j] * up
                ifault = 2
            end
        end
    @label endofproc
        #gc() 
        trace[7] = count
        res = qfval
        return
end
#qfc([1, 2, 3], [1, 2, 3], [1, 2, 3], 3, 3, 10.0, 9.0, 8, 7.0, 5, 4.0, 3.0, [1, 2, 3], false)