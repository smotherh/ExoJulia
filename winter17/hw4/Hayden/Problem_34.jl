#Practice with ttv_nplanet.jl
using LsqFit
using PyPlot

include("ttv_wrapper_fixp3.jl")
include("ttv_nplanet.jl")

TransitTime1 = vec(readdlm("ttv_planet1.txt"));
TransitTime2 = vec(readdlm("ttv_planet2.txt"));
TransitTimes = vcat(TransitTime1, TransitTime2);

Period1 = mean(TransitTime1[2:end]-TransitTime1[1:end-1])
Period2 = mean(TransitTime2[2:end]-TransitTime2[1:end-1])

X = linspace(1,length(TransitTimes), length(TransitTimes))

#mass_ratio, period, trans0, ecosw, esinw
Params1 = [1, Period1, 0, 0.01, 0.01];
Params2 = [.9, Period2, 0, 0.01, 0.01];
Params3 = [1, 0, 0.1, 0.1];
Params = vcat(Params1, Params2, Params3);

NumTrans1 = length(TransitTime1)
NumTrans2 = length(TransitTime2)
NumTrans = [NumTrans1,NumTrans2]

#TTVs = ttv_nplanet(2, 5, NumTrans, vcat(Params1, Params2))

Residual=1e300;
BestFit=0.;
BestPeriod3=0.;
fit=0;

for i = 4300:.1:4350
  global p3_cur=i;
  fit = curve_fit(ttv_wrapper_fixp3, X, TransitTimes, Params)
  if sum((fit.resid).^2)<Residual
    Residual = sum(fit.resid.^2);
    BestFit = fit.param;
    BestPeriod3=p3_cur;
  end
end

Planet1 = Dict("MassRatio"=>0., "Period"=>0., "Trans0"=>0., "ecosw"=>0., "esinw"=>0., "Name"=>"Venus");
Planet1["MassRatio"] = BestFit[1];
Planet1["Period"]    = BestFit[2];
Planet1["Trans0"]    = BestFit[3];
Planet1["ecosw"]     = BestFit[4];
Planet1["esinw"]     = BestFit[5];

Planet2 = Dict("MassRatio"=>0., "Period"=>0., "Trans0"=>0., "ecosw"=>0., "esinw"=>0., "Name"=>"Earth");
Planet2["MassRatio"] = BestFit[6];
Planet2["Period"]    = BestFit[7];
Planet2["Trans0"]    = BestFit[8];
Planet2["ecosw"]     = BestFit[9];
Planet2["esinw"]     = BestFit[10];

Planet3 = Dict("MassRatio"=>0., "Period"=>0., "Trans0"=>0., "ecosw"=>0., "esinw"=>0., "Name"=>"Jupiter");
Planet3["MassRatio"] = BestFit[11];
Planet3["Period"]    = BestPeriod3;
Planet3["Trans0"]    = BestFit[12];
Planet3["ecosw"]     = BestFit[13];
Planet3["esinw"]     = BestFit[14];

function PrintPlanet(Planet)
  print(string("Orbital parameters from TTV fitting for planet ", Planet["Name"],":\n"))
  print(string("  Period: ", Planet["Period"],"\n"))
  print(string("  Mass Ratio: ", Planet["MassRatio"],"\n"))
  print(string("  Transit Time t0: ", abs(Planet["Trans0"]%Planet["Period"]),"\n"))
  print(string("  ecos(w): ", Planet["ecosw"],"\n"))
  print(string("  esin(w): ", Planet["esinw"],"\n"))
end

PrintPlanet(Planet1)
PrintPlanet(Planet2)
PrintPlanet(Planet3)
