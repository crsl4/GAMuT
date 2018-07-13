function exp1(x)               
    if -0.5 * sum1 < -50.0
        return(0)
    else 
        return(exp(x))
    end
end

function errbd(u, cx)#cx is a pointer to a double 
    xconst = u * sigsq #sigsq is external var. double 
    sum1 = u * xconst
    u = 2.0 * u
    for j in (r-1):0 #r is external variable integer
        nj = n[j] #int     ##n(int), lb (double), and nc(double) are pointers
        lj = lb[j] #double
        ncj = nc[j] #double
        x = u * lj #double
        y = 1.0 - x #double 
        xconst = xconst + lj * (ncj / y + nj) / y
        sum1 = sum1 + ncj * (x/y)^2 + nj * x^2 /y +log(x) #don't know how to convert log
    end
    cx = xconst #????
    return(exp1(-0.5 * sum1))
end

function ctff(accx, upn) #upn is a pointer
    u2 = upn
    u1 = 0.0
    c1 = mean #extern. double
    if u2 > 0.0
        rb = 2.0 * lmax #extern. double
    else
        rb = 2.0 * lmin #extern. double
    end
    for u = u2 / (1.0 + u2 * rb)
        if errbd(u, c2) > accx #don't know how to handle this
            u1 = u2
            c1 = c2
            u2 = 2.0 * u2
        end
    end
    for u = (c1 - mean) / (c2 - mean)
        if u < 0.9
            u = (u1 + u2) / 2.0
            if (errbd(u / (1.0 + u * rb), xconst) > accx)
                u1 = u
                c1 = xconst
             else
                u2 = u
                c2 = xconst
             end
        end 
    end
    upn = u2 # don't know how to call and already called this
    return(c2)
end

function truncation(u, tausq)
    sum1 = 0.0
    prod2 = 0.0
    prod3 = 0.0
    s = 0.0
    sum2 = (sigsq +tausq) * u^2
    prod1 = 2.0 * sum2
    u = 2.0 * u
    for j in 0:r
        lj = lb[j]
        ncj = nc[j]
        nj = n[j]
        x = (u * lj)^2
        sum1 = sum1 + ncj * x / (1.0 + x)
        if x > 1.0
            prod2 = prod2 + nj * log(x)
            prod3 = prod3 + nj * log(x) #may be wrong log function 
            s = s + nj
        else
            prod1 = prod1 + nj * log(x) #may be wrong log
        end
    end
    sum1 = 0.5 *sum1 #this would stil be zero?!?!
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
    if err2 < err1 ||  err1 == err2
        x = 0.5 * sum2
    end
    err2 =  ( x  <=  y )  ? 1.0  : y / x
    return(err1 < err2) ? err1 : err2
end


function findu(utx, accx) ##no idea 
    divis = [2.0, 1.4, 1.2, 1.1]
    ut = utx
    u = ut / 4
    if truncation(u, 0.0) > accx
        ut = ut * 4.0
    else
        u = u /4
    end

end


function integrate(nterm, interv, tausq, mainx)
    
end

