function [gamma] = quartic_root (vrec, vtrans, Rs, author)

% Solve the correct root of quartic (gamma) for Martin-Neira (1993) and 
% Helm (2008) equations

xt=vtrans(1); % X coordinate of transmitter/satellite
yt=vtrans(2); % Y coordinate of transmitter/satellite

xr=vrec(1); % X coordinate of receiver
yr=vrec(2); % Y coordinate of receiver

%% Coefficients of quartic polynomial
c0 = ((xt.*yr)+(yt.*xr))-(Rs.*(yt+yr));
c1 = (-4.*((xt.*xr)-(yt.*yr)))+(2.*Rs.*(xt+xr));
c2 = (-6.*((xt.*yr)+(yt.*xr)));
c3 = (4.*((xt.*xr)-(yt.*yr)))+(2.*Rs.*(xt+xr));
c4 = ((xt.*yr)+(yt.*xr))+(Rs.*(yt+yr));

%% Starting estimate value 
cos_gamma_rt = dot(vtrans, vrec)./(norm(vtrans).*norm(vrec));
gamma_rt = acosd(cos_gamma_rt); % Angle between receiver and transmitter direction

gamma0 = (90-(gamma_rt./3));
t0 = tand(gamma0./2); % Starting estimate value 

%% Solution to gamma - Martín-Neira
% Solution to roots of quartic polynomial with roots.m function

if strcmp(author,'martin-neira')==1
    
    cs = [c4 c3 c2 c1 c0];
    ts = roots(cs); % Roots of quartic
    [~,ind] = min(abs(ts-t0));
    t = ts(ind);
    gamma = atand(t).*2; % Rotation angle phi, p. 333

end

%% Solution to gamma - Helm
% Solution to roots of quartic polynomial iterating gamma (without using all of roots)
% Iterative scheme based on modified Newthon-method - p.35-36

if strcmp(author,'helm')==1
    
    i = 1;
    i_max = 1000;
    tns(1) = tand(gamma0/2);
    gammas(1) = gamma0;
    % cond = 1e-8*(180/pi);
    cond = sqrt(eps());
    
    while (i <= i_max)
        
        tn = tns(i);
        fn = c4.*tn.^4 + c3.*tn.^3 + c2.*tn.^2 + c1.*tn + c0;
        f1n = 4.*c4.*tn.^3 + 3.*c3.*tn.^2 + 2.*c2.*tn + c1;
        f2n = 12.*c4.*tn.^2 + 6.*c3.*tn + 2.*c2;
        
        kn = abs(f2n./f1n);
        
        tns(i+1) = tn-(fn./f1n);
        gammas(i+1) = 2.*atand((tns(i+1)));
        
        if (abs((gammas(i+1))-(gammas(i))) < cond)
            break
        end
        
        if (f1n==0 && kn>1) == true
            break
        end
        i = i+1;
    end
    
    if (i > i_max)
        disp('Did not converged!');
    end
    
    gamma = gammas(i);
    
end