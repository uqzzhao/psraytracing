'''
@author:     Zhengguang Zhao
@copyright:  Copyright 2016-2019, Zhengguang Zhao.
@license:    MIT
@contact:    zg.zhao@outlook.com

'''
function psraytrace(vp, vs, zlayer, dg, sourcex, sourcey, sourcez, receiverx, receivery, receiverz)
   
    #   Input geometry
    max_x = maximum(sourcex)
    max_y = maximum(sourcey)
    topright_dis = ceil(sqrt(max_x * max_x + max_y * max_y)) + dg
    xmin = 0
    xmax = topright_dis


    #   Make strata layer
    xx = collect(xmin:dg:xmax); 


    #   Source-Receiver Groups
    #   Source
    zs = sourcez
    xs = sourcex - sourcex
    ns = length(xs)
    #     # Receiver
    zr = receiverz
    nr = length(zr)

    nray = ns * nr
    
    times = zeros(nray,1)
    tetas = zeros(nray,1)


    
    # Run Ray Tracing
       
    # Loop over for number of sources
    for i = 1: ns
        # Loop over for number of receiver
         
        for j = 1: nr
            
            xr[j] = sqrt((sourcex[i] - receiverx[j]) * (sourcex[i] - receiverx[j]) + (sourcey[i] - receivery[j]) * (sourcey[i] - receivery[j]))
            # Compare zs and zr and determine downgoing or upgoing

            if zs[i] > zr[j]

                # Upgoing path()
                u = findall(x->(x<zs[i]), zlayer)
                
                if isempty(u)
                    eup = length(zlayer)
                else
                    eup = u[end]
                end
                u = findall(x->(x>zr[j]), zlayer)
                if isempty(u)
                    sup = length(zlayer)
                else
                    sup = u[1]
                end
                zu = [zr[j];zlayer[sup:eup];zs[i]]
                nu = length(zu)
                zn = reverse(zu)

                # Upgoing elastic parameter
                vpu = vp[sup-1:eup]
                vsu = vs[sup-1:eup]

                # Combine model elastic parameter
                vpp = reverse(vpu)
                vps = reverse(vsu)

                # Start Raytracing [P-P, S-S, or P-S mode]
                ops = 1; # ops=1 for PP mode; ops=2 for PS mode
                
                xh,zh,vh,pp,teta,time = shooting(vpp,vps,zn,xx,xs[i],xr[j],ops)


            elseif zs[i] == zr[j]

                # Horizontal path()
                h = findall(x->(x<zs[i]), zlayer)
                
                if isempty(h)
                    hor = 1
                else
                    hor = h[end]
                end

                zhor = [zs[i] zr[j]]
                nu = length(zhor)
                zn = zhor

                # Upgoing elastic parameter
                vph = vp[hor]
                vsh = vs[hor]

                # Combine model elastic parameter
                vpp = vph
                vps = vsh

                # Start Raytracing [P-P, S-S, or P-S mode]
                ops = 1; # ops=1 for PP mode; ops=2 for PS mode
                xh,zh,vh,pp,teta,time = directshooting(vpp,vps,zn,xx,xs[i],xr[j],ops)



            else()
                # Downgoing path()
                d = findall(x->(x>zs[i]), zlayer)
                
                if isempty(d)
                    sdown = length(zlayer)
                else
                    sdown = d[1]
                end
                d = find(zlayer < zr[j])
                if isempty(d)
                    edown = length(zlayer)
                else
                    edown = d[end]
                end
                zd = [zs[i]; zlayer[sdown:edown]; zr[j]]
                nd = length(zd)
                zn = zd

                # Downgoing elastic parameter
                vpd = vp[sdown-1:edown]
                vsd = vs[sdown-1:edown]

                # Combine model elastic parameter
                vpp = vpd
                vps = vsd

                # Start Raytracing [P-P, S-S, or P-S mode]
                ops = 1; # ops=1 for PP mode; ops=2 for PS mode
                xh,zh,vh,pp,teta,time = shooting(vpp,vps,zn,xx,xs[i],xr[j],ops)


            end
            
            # Store traveltimes and incidence angles
            times[j,i] = time
            tetas[j,i] = teta[end]

            # Plot Ray: @TODO add return rays
            dis = 0 
            if dis == 1
                L = sqrt((sourcex[i] - receiverx[j]) * (sourcex[i] - receiverx[j]) + (sourcey[i] - receivery[j]) * (sourcey[i] - receivery[j]))
                X = sourcex[i] - receiverx[j]
                Y = sourcey[i] - receivery[j]

                if X <= 0
                    dx = sourcex[i] + xh/L * abs(X)
                else()
                    dx = sourcex[i] - xh/L * abs(X)
                end

                if Y <= 0
                    dy = sourcey[i] + xh/L * abs(Y)
                else()
                    dy = sourcey[i] - xh/L * abs(Y)
                end

            end

        end
    
    end
    
    times, tetas    
      
end


function shooting(vpp,vps,zn,xx,xs,xr,ops)

    # Some Constants
    itermax = 50
    offset = abs.(xs-xr)
    
    xc = 10.0
     
    xh = vcat([xs],xx)    
    time = 0.0
    teta = 0.0
    

    # Determine Option
    if ops == 1
        vh = vpp
    elseif ops == 2
        vh = vps
    end
    # Initial guess of the depth & time
    zh = zn #.- 100000*eps()
   
    t = Inf * ones(length(offset),1)
    p = Inf * ones(length(offset),1)

    # Start Raytracing
    # Trial shooting()
    pmax = 1/minimum(vh)
    pp = collect(range(0,1/maximum(vh),length = length(xx)))
    pp = reshape(pp, (1, length(pp)))
    
    
    sln = vh[1:length(zh)-1]*pp
    vel = vh[1:length(zh)-1]*ones(1,length(pp))    
    dz =  abs.(diff(zh))*ones(1,length(pp))
    
    if size(sln,1)>1
        xn = sum((dz.*sln)./sqrt.(1 .- sln.^2), dims = 1)
        tt = sum(dz./(vel.*sqrt.(1 .- sln.^2)), dims = 1)
        
    else()
        xn = (dz.*sln)./sqrt.(1 .- sln.^2)
        tt = dz./(vel.*sqrt.(1 .- sln.^2))
    end
    
    xmax = maximum(xn)   
    

    # Bisection Method
    # Start Bisection Method
    for k=1:length(offset)
        # Analyze the radius of target
        n = length(xn)
        
        xa = xn[1:n-1]
        xb = xn[2:n]        
        
        opt1 = (xa .<= offset[k]) .& (xb .> offset[k])        
        opt2 = (xa .>= offset[k]) .& (xb .< offset[k])
              
        opts = opt1 .| opt2          
        ind = findall(x->x==true, opts)
        
        if isempty(ind)
            if offset[k] >= xmax
                a = n
                b = []
            else
                a = []
                b = 1
            end
        else
            a = ind
            b = ind .+ 1
        end

        x1 = xn[a][1]
        x2 = xn[b][1]
        t1 = tt[a][1]
        t2 = tt[b][1]
        p1 = pp[a][1]
        p2 = pp[b][1]
        iter = 0
        err= ((b-a)/2)[1]
        
              
        # Minimize the error() & intersect the reflector
        while (iter < itermax) & (abs(err) < 1)
            iter = iter + 1
            println("iter： $iter")
            xt1 = abs.(offset[k] - x1)
            xt2 = abs.(offset[k] - x2)
            
            if (xt1 < xc) & (xt1 <= xt2)
                # Linear interpolation                
                t[k] = t1 + (offset[k] - x1)*(t2-t1)/(x2-x1)
                p[k] = p1 + (offset[k] - x1)*(p2-p1)/(x2-x1)
            elseif (xt2 < xc) & (xt2<=xt1)
                # Linear interpolation
                t[k] = t2 + (offset[k] - x2)*(t1-t2)/(x1-x2)
                p[k] = p2 + (offset[k] - x2)*(p1-p2)/(x1-x2)
            end
            # Set new ray parameter
            if isempty(a)
                p2 = p1
                p1 = 0
            elseif isempty(b)
                p1 = p2
                p2 = pmax
            end
            
            pnew = collect(range(min(p1,p2)[1],max(p1,p2)[1],length = 3))            
            pnew = reshape(pnew, (1, length(pnew)))            
            
            # Do shooting by new ray parameter            
            sln = vh[1:length(zh)-1]*pnew[2]            
            vel = vh[1:length(zh)-1]*ones(1,length(pnew[2]))
            dz = abs.(diff(zh))*ones(1,length(pnew[2]))
            if size(sln,1)>1
                xtemp = sum((dz.*sln)./sqrt.(1 .- sln.^2), dims = 1)
                ttemp = sum(dz./(vel.*sqrt.(1 .- sln.^2)), dims = 1)
            else
                xtemp = (dz.*sln)./sqrt.(1 .- sln.^2)
                ttemp = dz./(vel.*sqrt.(1 .- sln.^2))
            end
            xnew = [x1 xtemp x2]
            tnew = [t1 ttemp t2]
            xmax = maximum(xnew)
            # Analyze the radius of target
            n = length(xnew)
            xa = xnew[1:n-1]
            xb = xnew[2:n]
            opt1 = (xa .<= offset[k]) .& (xb .> offset[k])
            opt2 = (xa .>= offset[k]) .& (xb .< offset[k])
            
            opts = opt1 .| opt2          
            ind = findall(x->x==true, opts)
            
            a = ind
            b = ind .+ 1
            x1 = xnew[a][1]
            x2 = xnew[b][1]
            t1 = tnew[a][1]
            t2 = tnew[b][1]
            p1 = pnew[a][1]
            p2 = pnew[b][1]
            err = ((b - a)/2)[1]
            # Declare ray parameter
            if xr > xs
                pp = p
            else
                pp = -p
            end
            
            # Compute travel time & angle()            
            dx = real(complex.(pp.*vh.*dz)./sqrt.(complex.(1 .- pp.*pp.*vh.*vh)))
            xx = xs .+ cumsum(dx, dims = 1)                
            xh = vcat([xs],xx)  
            dz = real(dx.*sqrt.(complex.(1 .- pp.*pp.*vh.*vh))./(complex.(pp.*vh)))
            dt = dz./(vh.*sqrt.(complex.(1 .- pp.*pp.*vh.*vh)))

                       
            tt = cumsum(dt, dims =1)
            time = real(tt[end])
            println("time: $time")
            teta = real(asin.(complex.(pp.*vh)))          
         

        end # End the while loop
        
    end # End the offset loop
    println("zh: $zh")
    xh,zh,vh,pp,teta,time
    
end


function directshooting(vpp,vss,zn,xx,xs,xr,ops)

    # Horizontal path()

    if xs < xr
        xh = [xs xr]
    else
        xh = [xr xs]
    end

    zh = zn
    vh = vpp
    teta = 0.0

    if ops == 1
        pp = 1/vpp
        time = abs(xs-xr)/vpp
    else
        pp = 1/vss
        time = abs(xs-xr)/vss
    end

    xh,zh,vh,pp,teta,time
    
end # End the offset loop
