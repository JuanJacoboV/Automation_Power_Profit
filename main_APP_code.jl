using   CairoMakie,LaTeXStrings,NLboxsolve,Distributions,HPFilter,XLSX,ReadStat,DataFrames,CSV

##### This code reproduces the results in Section IV in the main text with the NAIRU.
## Calibrated (static) parameters
α = 1.64 # labor augmenting parameter
σ=0.6    # elasticity of substitution
Ak=0.02  # capital auugmenting parameter
A= 1     # Hicks neutral parameters
ρm=(1+0.26)^(1/12)-1   # this is obtained from Andreoni and Sprenger (2012)
gm = (1+0.02)^(1/12)-1   # this is the average of labor productivity growth
δm = (1+0.068)^(1/12)-1  # this is the average of the depreciation rate 
# Beveridge curve
λ0 = 0.025  # Separation rate
ṁ = 0.0000  # Equilibrium value
Ṁ = gm/α  # derived from the steady-state growth equation. 
ι = 1.25      #  Same as Petrosky-Nadeu (2018)
ρ = ρm    # subjective discount rate
γf = 0.45   #  Response time of capitalists
ξ = 8.2    # calibrated to be about 2 times the average productivity of labor
δ = δm  # obtained from Barki 
φ = 0.3   # Relative price of capital

# Code for solving the steady-state (with Tw as an endogenous variable)
function sys2!(FF,x,m,α,σ,λ0,ι,ρ,ṁ,Ṁ,γf,b_y,ξ,δ,L,Pr,A,Ak,φ)
	w=x[1]; wu=x[2]; wn=x[3];  k=x[4]; ku=x[5]; kn=x[6];
	μ=x[7]; θ = x[8]; θu = x[9]; θn = x[10]; fθ = x[11]; qθ = x[12];
	fθu = x[13]; qθu = x[14]; fθn= x[15]; qθn = x[16]; Tw = x[17];
	Lu=x[18]; Ln=x[19]; Γna = x[20]; Γnb = x[21]; Ψna = x[22]; 
	Ψnb=x[23] ; Ψn = x[24]; y = x[25]; ∂yL = x[26]; ∂yk = x[27];
	yu = x[28]; ∂yLu=x[29]; yn =x[30]; ∂yLn=x[31]; ∂ykn=x[32];∂yku=x[33];
	μu = x[34]; μn = x[35]; Ψu = x[36]; UAL = x[37]; g = x[38]; λ = x[39]; b=x[40]
	FF[1] = qθ - (1+θ^(ι))^(-1/ι)        #  probability of filling a vacancy
	FF[2] = fθ - qθ*θ                    # job finding probability
	FF[3] = L - fθ/(λ+fθ)                # steady-state employment
	FF[4] = qθu - (1+θu^(ι))^(-1/ι)      # probability of filling a vacancy in collective bargaining
	FF[5] = fθu - qθu*θu                 # job finding probability in collective bargaining
	FF[6] = Lu - fθu/(λ+fθu)             # steady-state employment in collective bargaining
	FF[7] = qθn - (1+θn^(ι))^(-1/ι)      # probability of filling a vacancy in individual bargaining
	FF[8] = fθn - qθn*θn                 # job finding probability in individual bargaining
	FF[9] = UAL - (1-exp(α*(σ-1)*(Ṁ-ṁ))*(exp(α*(σ-1)*(m+ṁ))-1)/(exp(α*(σ-1)*m)-1)) # marginal technological unemployment
	FF[10] = λ - λ0 - UAL                # separation rate
	FF[11] = g - α*Ṁ                     # growth rate
	FF[12] = Ln- fθn/(λ+fθn)             # steady-state employment in individual bargaining
	FF[13] = Γna - γf/(1+γf)             # Intrinsic bargaining power (case a)
	FF[14] = Γnb-γf*(1-qθn)/(1+γf+qθn*(1-γf)) # Intrinsic bargaining power (case b)
	FF[15] = Ψna - Γna*(ρ-g+λ+fθn)/(ρ-g+λ+	Γna*fθn)  # Actual bargaining power (case a)
	FF[16] = Ψnb - Γnb*(ρ-g+λ+fθn)/(ρ-g+λ+	Γnb*fθn)  # Actual bargaining power (case b)
	FF[17] = Ψn - ((Tw*Γnb +θn*Γna)/(θn+Tw))*(ρ-g+λ+fθn)/(ρ-g+λ+((Tw*Γnb +θn*Γna)/(θn+Tw))*fθn) # Actual bargaining power (individual protocol)
	FF[18] = Ψu - Γna*(ρ-g+λ+fθu)/(ρ-g+λ+	Γna*fθu)  # Actual bargaining power (collective protocol)
	FF[19] = y -A*((1-m)^(1/σ) *(Ak*k)^((σ-1)/σ) + ((exp(α*(σ-1)*m)-1)/(α*(σ-1)))^(1/σ))^(σ/(σ-1)) # production function
	FF[20] = ∂yk - (Ak*A)^((σ-1)/σ)*(y/k)^(1/σ) * (1-m)^(1/σ)        # Marginal prod. capital 
	FF[21] = ∂yL - (y - k*∂yk)                                   # Marginal prod. labor 
	FF[22] = yu - A*((1-m)^(1/σ) *(Ak*ku)^((σ-1)/σ) + ((exp(α*(σ-1)*m)-1)/(α*(σ-1)))^(1/σ))^(σ/(σ-1))  # production function collective bargaining
	FF[23] = ∂yku - (A*Ak)^((σ-1)/σ)*(yu/ku)^(1/σ) * (1-m)^(1/σ)          # Marginal prod. capital with collective bargaining
	FF[24] = ∂yLu - (yu - ku*∂yku)                                        #      Marginal prod. labor with collective bargaining
	FF[25] = yn - A*((1-m)^(1/σ) *(Ak*kn)^((σ-1)/σ) + ((exp(α*(σ-1)*m)-1)/(α*(σ-1)))^(1/σ))^(σ/(σ-1)) # production function individual bargaining
	FF[26] = ∂ykn - (A*Ak)^((σ-1)/σ)*(yn/kn)^(1/σ) * (1-m)^(1/σ)                # Marginal prod. capital (equation (A1)) with individual bargaining
	FF[27] = ∂yLn -(yn - kn*∂ykn)                                        # Marginal prod. labor (equation (A1)) with individual bargaining
	FF[28] = wn -( b + Ψn*(∂yLn-b))                                      # Hypothetical wage individual bargaining 
	FF[29] = wn -( ∂yLn- ξ*(ρ-g+λ)/qθn)                                  # Hypothetical labor demand with individual bargaining  (these equations solve θn,wn)
	FF[30] = wu -( b + Ψu*(∂yLu-b + (ρ-g+λ)/(ρ-g) *(yu-∂yLu)))           # Hypothetical wage collective bargaining 
	FF[31] = wu -( ∂yLu- ξ*(ρ-g+λ)/qθu )                                 # Hypothetical labor demand with collective bargaining  (these equations solve θu,wu)
	FF[32] = w -(Pr*wu + (1-Pr)*wn)                                      # Aggregate wage   
	FF[33] = w -(∂yL- ξ*(ρ-g+λ)/qθ)                                      # Labor demand equation   
	FF[34] = μ -((∂yL/w)-1)                                              # markup
	FF[35] = ∂yk - δ*(1+μ)/φ                                             # Capital market equilibrium
	FF[36] = μn -((∂yLn/wn)-1)                                           # markup (individual bargaining)
	FF[37] = μu - ((∂yLu/wu)-1)                                          # markup (collective bargaining)
	FF[38] = ∂yku - δ*(1+μu)/φ                                           # Capital market equilibrium (collective bargaining)
	FF[39] = ∂ykn - δ*(1+μn)/φ                                           # Capital market equilibrium (individual bargaining)
	FF[40] = b - b_y*y                                                   # Steady-state unemployment benefits
end


## Time-Varying Parameters
# Download the necessary data
data_calibration = XLSX.readxlsx("/Users/juanjacobo/Documents/Documents/metroeconomica_2025/Data/data_in_calibrations.xlsx");
data_annual=data_calibration["annual!A1:AW83"];
data_quarterly=data_calibration["quarterly!A1:N290"]
data_monthly=data_calibration["monthly!A1:G962"]
## the data is annual from 1948 to 2008
time_annual = convert(Array{Float64,1},data_annual[12:end-11,1])
g_penn = convert(Array{Float64,1},data_annual[12:end-11,5])  # growth rate penn world data
unions = convert(Array{Float64,1},data_annual[12:end-11,2])  # union data
poor_welfare = convert(Array{Float64,1},data_annual[26:end-11,45])  # 1962-2010
poor_welfare_PTANF = convert(Array{Float64,1},data_annual[26:end-11,46])  # 1962-2010
poor_welfare_PTANFMEDICAID = convert(Array{Float64,1},data_annual[26:end-11,47])  # 1962-2010
poor_welfare_ALL_LTAXCREDITS = convert(Array{Float64,1},data_annual[26:end-11,48])  # 1962-2010
automation_mann = convert(Array{Float64,1},data_annual[38:72,49])   # 1974-2008

## the data is quarterly from 1948 to 2008
time_quarterly = collect(range(1948,length= 245, stop=2009)) 
b_chodorow_q = convert(Array{Float64,1},data_quarterly[6:end-40,4]) # Opportunity costs from Chodorow-Reich, & Karabarbounis
NRU_q= convert(Array{Float64,1},data_quarterly[10:end-40,5])        # Noncyclical rate of unemployment
YL=convert(Array{Float64,1},data_quarterly[6:end-40,2])             # GDP per capita
# Compute moving averages
b_chodorow =zeros(61)
NRU = [5.25;zeros(60)]
yl = zeros(61)
for i=1:61
	b_chodorow[i]=mean(b_chodorow_q[1+(i-1)*4:i*4])
	yl[i]=mean(YL[1+(i-1)*4:i*4])
end
for i=1:60
	NRU[i+1]=mean(NRU_q[1+(i-1)*4:i*4]) # Natural rate of unemployment
end
## the data is monthly from 1948 to 2009
time_monthly = collect(range(1948,length= 733, stop=2009)) 
un_civil_m = convert(Array{Float64,1},data_monthly[110:842,5]) # civil unemployment
vacancies_m = convert(Array{Float64,1},data_monthly[110:842,7]) # monthly vacancies
un_c,vac=zeros(61),zeros(61),zeros(61),zeros(61)
for i=1:61
	un_c[i]=mean(un_civil_m[1+(i-1)*12:i*12])
	vac[i]=mean(vacancies_m[1+(i-1)*12:i*12])
end
θs = vac./un_c   # Labor market tightness
bcho_HP = HP(b_chodorow,6)  # unemployment benefts (b)
unions_HP = HP(unions,6)  # Union data
G_hp = HP(g_penn,6)  # Penn world Table growth data
g_ms = (G_hp.+1).^(1/12).-1 # monthly 
Ṁs = 1.0*g_ms./α;  # Growth in the creation of new tasks


#bea_bls_data is the data from Eldridge et al. (2020)
data = XLSX.readxlsx("/Users/juanjacobo/Documents/Documents/metroeconomica_2025/Data/bea_bls_data.xlsx")
data_depre = XLSX.readxlsx("/Users/juanjacobo/Documents/Documents/metroeconomica_2025/Data/current_dep_assets.xlsx")
curr_c_cap = XLSX.readxlsx("/Users/juanjacobo/Documents/Documents/metroeconomica_2025/Data/cost_assets.xlsx")
data_1947_1963=data["1947-1963!A2:W750"]
data_1963_2016=data["1963-2016!A2:W3404"]
full_sample=[data_1947_1963;data_1963_2016] # merges the time periods
data_D=data_depre["Sheet2!A1:U72"] # depreciation data
cost_assets=curr_c_cap["Sheet2!A1:U72"] # value of assets
# The following computes the empirical variables in the text using the experimental BEA-BLS data
function economy_lab_shares(full_sample,data_D,cost_assets)
	VA_t=zeros(17);  		VA_t2=zeros(54);  		VA_t3=zeros(71);  	
	la_c_t=zeros(17);		la_c_t2=zeros(54);  		la_c_t3=zeros(71);  	
	la_nc_t=zeros(17); 	la_nc_t2=zeros(54);  		la_nc_t3=zeros(71);  	
	MA_t=zeros(17);     	MA_t2=zeros(54);     MA_t3=zeros(71)
	KA_t=zeros(17);     	KA_t2=zeros(54);     KA_t3=zeros(71)
	rVA_t=zeros(17);  		rVA_t2=zeros(54);  		rVA_t3=zeros(71);  	
	L_t = zeros(17);  		L_t2=zeros(54);  		L_t3=zeros(71);  	
	h_t=zeros(17);  		h_t2=zeros(54);  		h_t3=zeros(71);  	
	depret=sum(data_D[2:end-1,i]./1 for i ∈ [5,6,7,8,9,10,11,16,19,20,21]).*1000
	Kts=sum(cost_assets[2:end-1,i]./1 for i ∈ [5,6,7,8,9,10,11,16,19,20,21]).*1000
	for i =[6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,59,60,61]
		agr_VA=findall(x -> x ==i, full_sample[:,2])
		VA=full_sample[agr_VA,4].-full_sample[agr_VA,5]
		rVA=full_sample[agr_VA,13].-full_sample[agr_VA,14]
		LL = (full_sample[agr_VA,20].+full_sample[agr_VA,21]).*full_sample[agr_VA,22]
		VA_t3=VA_t3.+VA
		rVA_t3=rVA_t3.+rVA
		MA_t3=MA_t3.+full_sample[agr_VA,5]
		la_c=full_sample[agr_VA,11]
		la_c_t3=la_c_t3+la_c
		la_nc=full_sample[agr_VA,12]
		h_nc=full_sample[agr_VA,22]
		la_nc_t3=la_nc_t3+la_nc
		h_t3=h_nc+h_t3
		L_t3 = L_t3 + LL
	end
	for i=[2936,3740,5152,5758]
		agr_VA=findall(x -> x ==i, full_sample[:,2])
		VA=full_sample[agr_VA,4].-full_sample[agr_VA,5]
		rVA=full_sample[agr_VA,13].-full_sample[agr_VA,14]
		VA_t=VA_t.+VA
		rVA_t=rVA_t.+rVA
		LL = (full_sample[agr_VA,20].+full_sample[agr_VA,21]).*full_sample[agr_VA,22]
		MA_t=MA_t.+full_sample[agr_VA,5]
		la_c=full_sample[agr_VA,11]
		la_c_t=la_c_t+la_c
		la_nc=full_sample[agr_VA,12]
		la_nc_t=la_nc_t+la_nc
		h_nc=full_sample[agr_VA,22]
		h_t=h_nc+h_t
		L_t = L_t + LL
	end
	for i=[29,30,31,32,33,34,35,36,37,38,39,40,51,52,57,58]
		agr_VA=findall(x -> x ==i, full_sample[:,2])
		VA=full_sample[agr_VA,4].-full_sample[agr_VA,5]
		rVA=full_sample[agr_VA,13].-full_sample[agr_VA,14]
		LL = (full_sample[agr_VA,20].+full_sample[agr_VA,21]).*full_sample[agr_VA,22]
		VA_t2=VA_t2.+VA
		rVA_t2=rVA_t2.+rVA
		MA_t2=MA_t2.+full_sample[agr_VA,5]
		la_c=full_sample[agr_VA,11]
		la_c_t2=la_c_t2+la_c
		la_nc=full_sample[agr_VA,12]
		la_nc_t2=la_nc_t2+la_nc
		h_nc=full_sample[agr_VA,22]
		h_t2=h_nc+h_t2
		L_t2 = L_t2 + LL
	end
	la_c=[la_c_t;la_c_t2].+la_c_t3; 
	la_nc=[la_nc_t;la_nc_t2].+la_nc_t3; 
	L_total =[L_t;L_t2].+L_t3
	h = [h_t;h_t2].+h_t3; 
	VA=[VA_t;VA_t2].+VA_t3
	rVA= [rVA_t;rVA_t2].+rVA_t3
	YL = VA./L_total
	PL = rVA./(h.*(la_c.+la_nc))
	MA=[MA_t;MA_t2].+MA_t3
	KA= [KA_t;KA_t2].+KA_t3
	Cap=la_nc[2:end].+la_c[2:end].+depret.+MA[2:end]
	Cap2 = la_nc[2:end].+la_c[2:end].+depret
	Prof=(VA[2:end].-depret.-la_nc[2:end].-la_c[2:end])
	VA_N=VA[2:end].-depret
return [(la_nc)./VA (la_nc.+la_c)./VA [0;Prof./Cap2] [0;(la_nc[2:end].+la_c[2:end])./VA_N] [0;Prof./(Kts.+depret)] [0;(Kts.+depret)./VA[2:end]] [0;(Prof.+depret.+la_nc[2:end].+la_c[2:end])./(VA[2:end])] [0;depret./(Kts.+depret)]]
end
economy_LS=economy_lab_shares(full_sample,data_D,cost_assets)
LS_NCP =economy_LS[:,1]         # labor share (unskilled labor)
LS_C = economy_LS[:,2].-LS_NCP  # labor share (skilled labor)
w_share =  economy_LS[1:end,2]  # labor share 
μ_data =  economy_LS[2:end,3]   # rate of return
r_data =  economy_LS[2:end,5]    # rate of profit
KY = economy_LS[2:end,6] ;        # capital-output ratio
dep_cap=economy_LS[2:end,end]; # depreciation  of capital

KYm=KY.*12  # remember that the calibration is based on monthly data
dep_capm = (dep_cap.+1).^(1/12).-1
m_star =-((mean(dep_capm).*(μ_data.+1)).^(σ)* (Ak*φ)^(1-σ).*(KYm)).+1 
mm=-m_star.+1;


#automation_hemous.csv is obtained from hemous_2025.jl
merged_ctry=CSV.read("/Users/juanjacobo/Documents/Documents/metroeconomica_2025/Code/automation_hemous.csv", DataFrame)

## Automation from Hemous 2025. Auto95 US
Mech=merged_ctry.Auto95_USA # data  1963:2015
Xmech=[ones(46) Mech[1:end-7]]
M̂1 = Xmech*(inv(Xmech'Xmech)*(Xmech'mm[16:end-9]))
mM̂2=-M̂1.+1;

### This code sets the initial values. 
function init_v(m,α,σ,λ0,ι,ρ,ṁ,Ṁ,γf,b_y,ξ,δ,Tw,Pr,A,Ak)
	k0=10.95; ku0=11.3; kn0=10.65; θ0=0.7; θu0=0.016; θn0=0.72; g=α*Ṁ ;   UAL =(1-exp(α*(σ-1)*(Ṁ-ṁ))*(exp(α*(σ-1)*(m+ṁ))-1)/(exp(α*(σ-1)*m)-1));   λ = λ0 + UAL; 
	qθ0=(1+θ0^(ι))^(-1/ι); fθ0=qθ0*θ0 ;   qθu0= (1+θu0^(ι))^(-1/ι); fθu0 = qθu0*θu0; qθn0=(1+θn0^(ι))^(-1/ι);  fθn0=qθn0*θn0;  L0=fθ0/(fθ0+λ); Lu0=fθu0/(fθu0+λ);
	Ln0=fθn0/(fθn0+λ); Γna0 = γf/(1+γf);  Γnb0 = γf*(1-qθn0)/(1+γf+qθn0*(1-γf)); 
	Ψna0=Γna0*(ρ-g+λ+fθn0)/(ρ-g+λ+	Γna0*fθn0); 
	Ψnb0=Γnb0*(ρ-g+λ+fθn0)/(ρ-g+λ+	Γnb0*fθn0); 
	Ψ0 = ((Γna0*θ0+Γnb0*Tw)/(Tw+θ0))*(ρ-g+λ+fθn0)/(ρ-g+λ+	((Γna0*θ0+Γnb0*Tw)/(Tw+θ0))*fθn0); 
	Ψu0=Γna0*(ρ-g+λ+fθu0)/(ρ-g+λ+	Γna0*fθu0);
	y0 = A*((1-m)^(1/σ) *(Ak*k0)^((σ-1)/σ) + ((exp(α*(σ-1)*m)-1)/(α*(σ-1)))^(1/σ))^(σ/(σ-1))
	b = b_y*y0
	∂yk0= (Ak*A)^((σ-1)/σ)*(y0/k0)^(1/σ) * (1-m)^(1/σ)
	∂yL0 = y0 - k0*∂yk0
	yu0=A*((1-m)^(1/σ) *(Ak*ku0)^((σ-1)/σ) + ((exp(α*(σ-1)*m)-1)/(α*(σ-1)))^(1/σ))^(σ/(σ-1))
	yn0=A*((1-m)^(1/σ) *(Ak*kn0)^((σ-1)/σ) + ((exp(α*(σ-1)*m)-1)/(α*(σ-1)))^(1/σ))^(σ/(σ-1))
	∂yku0 = (Ak*A)^((σ-1)/σ)*(yu0/ku0)^(1/σ) * (1-m)^(1/σ)
	∂ykn0 = (Ak*A)^((σ-1)/σ)*(yn0/kn0)^(1/σ) * (1-m)^(1/σ)
	∂yLu0 = yu0 - ku0*∂yku0
	∂yLn0 = yn0 - kn0*∂ykn0
	wu0=∂yLu0- ξ*(ρ-g+λ)/qθu0
	wn0=∂yLn0- ξ*(ρ-g+λ)/qθn0
	w0=∂yL0- ξ*(ρ-g+λ)/qθ0
	μu0=(∂yLu0/wu0)-1
	μn0=(∂yLn0/wn0)-1
	μ0=(∂yL0/w0)-1
	return [w0,wu0,wn0,k0,ku0,kn0,μ0,θ0,θu0,θn0,fθ0,qθ0,fθu0,
	qθu0,fθn0,qθn0,L0,Lu0,Ln0,Γna0,Γnb0,Ψna0,Ψnb0,Ψ0,y0,∂yL0,
	∂yk0,yu0,∂yLu0,yn0,∂yLn0,∂ykn0,∂yku0,μu0,μn0, Ψu0,UAL,g,λ,b]
	end

	lbt=[0.15,0.26,0.05,0.0,0.0,0.0,0,0,0.0005,0,0,0,0,0,0,0,0.35,
	0.0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.0,0,-1,-0.2,0,0] # lower bounds
	# upper bounds 
	ubt=[50,50,50,50,50,50,50,50,50,50,
	1,1,1,1.0001,1,1,6.5,1,1,1,1,1,1,1,
	10,10,10,10,10,10,10,10,10,10,10.0,1,0.1,0.1,0.1,1]

## Calibrate steady-state values assuming you know the NRU
function inverse_calibration(Pr,m,b_y,L,Ṁ,γf,λ0)
	x0=init_v(m,α,σ,λ0,ι,ρ,ṁ,Ṁ,γf,b_y,ξ,δ,1.5,Pr,A,Ak)
	x0[17]=1.75  # this is just an initial guess
	rlb = nlboxsolve((FF,x) ->sys2!(FF,x,m,α,σ,λ0,ι,ρ,ṁ,Ṁ,γf,b_y,ξ,δ,L,Pr,A,Ak,φ),
  x0,lbt,ubt,xtol=1e-4,ftol=1e-4)
  sol =rlb.zero
 # aggregate wage; capital; rate of exploitation; labor-market tightness;
 # rel. mobility of labor; output; worker power
 w=sol[1]; k=sol[4];μ = sol[7];θ= sol[8];Tw =sol[17];y =sol[25];Ψn =sol[24];b=sol[40]
	return [w;k;μ;θ;Tw;y;Ψn;b]
end
L_NRU=0.01*(-NRU.+100) # natural rate of unemployment
results_5H=zeros(8,46)
## Hypothesis: Changes in labor institutions  and technical change with NRU given
for i=1:46
 results_5H[:,i]=inverse_calibration(unions_HP[i+15],mM̂2[i],bcho_HP[i+15],
 L_NRU[i+15],Ṁs[i+15],γf,λ0) 
end
# Results
Ωs5H=results_5H[1,:]./results_5H[6,:] # Labor share
kys5H= results_5H[2,:]./results_5H[6,:] * (1/φ) # capital-output
Vs5H = results_5H[4,:].*NRU[16:end]./100 # vacancy rate
tsw5H = results_5H[5,:]                # Threat of competition
μs5H = results_5H[3,:]                 # markup
r5H = (μs5H./(μs5H.+1))./kys5H         # rate of profit
Ψn5H = results_5H[7,:]                 # Worker power
θ5sH = results_5H[4,:];                # labor market tightness


# Unlike the previous hypotheses, the rate of unemployment is determined by the model and Tw is an exogenous parameter  
function sys2d!(FF,x,m,α,σ,λ0,ι,ρ,ṁ,Ṁ,γf,b_y,ξ,δ,Tw,Pr,A,Ak,φ)
	w=x[1]; wu=x[2]; wn=x[3];  k=x[4]; ku=x[5]; kn=x[6];
	μ=x[7]; θ = x[8]; θu = x[9]; θn = x[10]; fθ = x[11]; qθ = x[12];
	fθu = x[13]; qθu = x[14]; fθn= x[15]; qθn = x[16]; L = x[17];
	Lu=x[18]; Ln=x[19]; Γna = x[20]; Γnb = x[21]; Ψna = x[22]; 
	Ψnb=x[23] ; Ψn = x[24]; y = x[25]; ∂yL = x[26]; ∂yk = x[27];
	yu = x[28]; ∂yLu=x[29]; yn =x[30]; ∂yLn=x[31]; ∂ykn=x[32];∂yku=x[33];
	μu = x[34]; μn = x[35]; Ψu = x[36]; UAL = x[37]; g = x[38]; λ = x[39]; b=x[40]
	FF[1] = qθ - (1+θ^(ι))^(-1/ι)
	FF[2] = fθ - qθ*θ
	FF[3] = L - fθ/(λ+fθ)
	FF[4] = qθu - (1+θu^(ι))^(-1/ι)
	FF[5] = fθu - qθu*θu
	FF[6] = Lu - fθu/(λ+fθu)
	FF[7] = qθn - (1+θn^(ι))^(-1/ι)
	FF[8] = fθn - qθn*θn
	FF[9] = UAL-(1-exp(α*(σ-1)*(Ṁ-ṁ))*(exp(α*(σ-1)*(m+ṁ))-1)/(exp(α*(σ-1)*m)-1))
	FF[10] = λ - λ0 - UAL
	FF[11] = g - α*Ṁ
	FF[12] = Ln- fθn/(λ+fθn)
	FF[13] = Γna - γf/(1+γf)
	FF[14] = Γnb-γf*(1-qθn)/(1+γf+qθn*(1-γf))
	FF[15] = Ψna - Γna*(ρ-g+λ+fθn)/(ρ-g+λ+	Γna*fθn)
	FF[16] = Ψnb - Γnb*(ρ-g+λ+fθn)/(ρ-g+λ+	Γnb*fθn)
	FF[17] = Ψn - ((Tw*Γnb +θn*Γna)/(θn+Tw))*(ρ-g+λ+fθn)/(ρ-g+λ+((Tw*Γnb +θn*Γna)/(θn+Tw))*fθn)
	FF[18] = Ψu - Γna*(ρ-g+λ+fθu)/(ρ-g+λ+	Γna*fθu)
	FF[19] = y -A*((1-m)^(1/σ) *(Ak*k)^((σ-1)/σ) + ((exp(α*(σ-1)*m)-1)/(α*(σ-1)))^(1/σ))^(σ/(σ-1))
	FF[20] = ∂yk - (Ak*A)^((σ-1)/σ)*(y/k)^(1/σ) * (1-m)^(1/σ)
	FF[21] = ∂yL - (y - k*∂yk)
	FF[22] = yu - A*((1-m)^(1/σ) *(Ak*ku)^((σ-1)/σ) + ((exp(α*(σ-1)*m)-1)/(α*(σ-1)))^(1/σ))^(σ/(σ-1))
	FF[23] = ∂yku - (A*Ak)^((σ-1)/σ)*(yu/ku)^(1/σ) * (1-m)^(1/σ)
	FF[24] = ∂yLu - (yu - ku*∂yku)
	FF[25] = yn - A*((1-m)^(1/σ) *(Ak*kn)^((σ-1)/σ) + ((exp(α*(σ-1)*m)-1)/(α*(σ-1)))^(1/σ))^(σ/(σ-1))
	FF[26] = ∂ykn - (A*Ak)^((σ-1)/σ)*(yn/kn)^(1/σ) * (1-m)^(1/σ)
	FF[27] = ∂yLn -(yn - kn*∂ykn)
	FF[28] = wn -( b + Ψn*(∂yLn-b))
	FF[29] = wn -( ∂yLn- ξ*(ρ-g+λ)/qθn)
	FF[30] = wu -( b + Ψu*(∂yLu-b + (ρ-g+λ)/(ρ-g) *(yu-∂yLu)))
	FF[31] = wu -( ∂yLu- ξ*(ρ-g+λ)/qθu )
	FF[32] = w -(Pr*wu + (1-Pr)*wn)
	FF[33] = w -(∂yL- ξ*(ρ-g+λ)/qθ)
	FF[34] = μ -((∂yL/w)-1)
	FF[35] = ∂yk - δ*(1+μ)/φ
	FF[36] = μn -((∂yLn/wn)-1)
	FF[37] = μu - ((∂yLu/wu)-1)
	FF[38] = ∂yku - δ*(1+μu)/φ
	FF[39] = ∂ykn - δ*(1+μn)/φ
	FF[40] = b - b_y*y
end
# lower bounds
lbd=[0.15,0.26,0.05,0.0,0.0,0.0,0,0,0.005,0,0,0,0,0,0,0,0,0.0,0,0,0,
0,0,0,0,0,0,0,0,0,0,0,0,0,0.0,0,-λ0,-0.2,0,0]
# upper bounds
ubd=[50,50,50,50,50,50,50,50,50,50,
1,1,1,1.0001,1,1,1,1,1,1,1,1,1,1,
10,10,10,10,10,10,10,10,10,10,10.0,1,0.1,0.1,0.1,1]

function inverse_calibrationd(Pr,m,b_y,Tw,Ṁ,γf,λ0,φ)
x0d=init_v(m,α,σ,λ0,ι,ρ,ṁ,Ṁ,γf,b_y,ξ,δ,Tw,Pr,A,Ak)
rlb = nlboxsolve((FF,x) ->sys2d!(FF,x,m,α,σ,λ0,ι,ρ,ṁ,Ṁ,γf,b_y,ξ,δ,Tw,Pr,A,Ak,φ),
x0d,lbd,ubd,xtol=1e-4,ftol=1e-4)
sol =rlb.zero
w=sol[1]; k=sol[4];μ = sol[7]; θ= sol[8]; L = sol[17]; y = sol[25]; Ψn = sol[24]; b=sol[40]
return [w;k;μ;θ;L;y;Ψn;b]
end
# Estimate equation (25)
X_tw=[ones(46) log.(((poor_welfare[1:end-1]))) log.(yl[15:end-1])]
β_tw=inv(X_tw'X_tw)*(X_tw'log.(tsw5H))
tsw_indiriect=X_tw*β_tw# sol[1].*cdf(Normal(0,1), X_tw*β̂_b)

with_theme(theme_latexfonts()) do       
    f = Figure(size = (500, 300))
    ax = Axis(f[1, 1],
        xlabel = "Period ",xgridstyle = :dash, ygridstyle = :dash,
        title = "Threat of competition among workers",
    )
    lines!(ax,time_annual[16:end],exp.(HP(tsw_indiriect,1000)), label = "Hypothesis")
    lines!(ax, time_annual[16:end],tsw5H, label = "Ideal")
    axislegend(; position = :lt, nbanks = 2, framecolor = (:grey, 0.95))
    f
    end

	#Automation alone: Task-displacing and labor-augmenting technical change
	results_1dH=zeros(8,46)
	for i=1:46
		results_1dH[:,i]=inverse_calibrationd(mean(unions_HP[16:end]),mM̂2[i],
		mean(bcho_HP[16:end]),mean(tsw5H),Ṁs[i+15],γf,λ0,φ)
	end
	# Automation and labor augmneting technical change
	ΩsdH=results_1dH[1,:]./results_1dH[6,:]
	kysdH= results_1dH[2,:]./results_1dH[6,:] * (1/φ)
	VsdH = results_1dH[4,:].*(-results_1dH[5,:].+1)
	UsdH = -results_1dH[5,:].+1
	μsdH = results_1dH[3,:]
	θsdH= results_1dH[4,:]
	rsdH = (12*μsdH./(μsdH.+1))./kysdH
	ΨndH= results_1dH[7,:];
	
	results_2dH=zeros(8,46)
	# Institutions and automation
	for i=1:46
		results_2dH[:,i]=inverse_calibrationd(unions[i+15],mM̂2[i],b_chodorow[i+15],
		exp.(tsw_indiriect)[i],Ṁs[i+15],γf,λ0,φ)
	end
	Ωsd2H=results_2dH[1,:]./results_2dH[6,:]
	kysd2H= results_2dH[2,:]./results_2dH[6,:] * (1/φ)
	Vsd2H = results_2dH[4,:].*(-results_2dH[5,:].+1)
	Usd2H = -results_2dH[5,:].+1
	μsd2H = results_2dH[3,:]
	θsd2H= results_2dH[4,:]
	rsd2H = (12*μsd2H./(μsd2H.+1))./kysd2H
	Ψnd2H= results_2dH[7,:];
	
	# Graphs
	
	with_theme(theme_latexfonts()) do       
		f = Figure(size = (900, 600))
		ax = Axis(f[1, 1],
			xlabel = " ",xgridstyle = :dash, ygridstyle = :dash,
			title = "Labor share",
		)
		lines!(ax, time_annual[16:end], ΩsdH, label="Automation")
		   lines!(ax, time_annual[16:end], HP(Ωsd2H,1000), label="Auto. and inst.")
			  lines!(ax, time_annual[16:end], w_share[16:end-10], label="Data")
			   axislegend(; position = :lb, nbanks = 1, framecolor = (:gray, 0.5))
		ax = Axis(f[1, 2],
		xlabel = " ",xgridstyle = :dash, ygridstyle = :dash,
		title = "Markup",
	)
	 lines!(ax, time_annual[16:end], μsdH)
		   lines!(ax, time_annual[16:end], HP(μsd2H,1000))
			  lines!(ax, time_annual[16:end], μ_data[17:end-8])
	  ax = Axis(f[2, 1],
	  xlabel = "Period ",xgridstyle = :dash, ygridstyle = :dash,
	  title = "Capital-output ratio",
	)
	lines!(ax, time_annual[16:end], kysdH./12)
		   lines!(ax, time_annual[16:end], kysd2H./12)
			  lines!(ax, time_annual[16:end],KY[16:end-9])
	ax = Axis(f[2, 2],
	xlabel = "Period",xgridstyle = :dash, ygridstyle = :dash,
	title = "Equilibrium Unemployment",
	)
	 lines!(ax, time_annual[16:end], UsdH)
		   lines!(ax, time_annual[16:end], HP(Usd2H,1000))
			  lines!(ax, time_annual[16:end], NRU[16:end]./100)
		f
		end
#

#CAGR
## Labor share
# 1963-1977
time_annual[16:30]
100*((HP(Ωsd2H,1000)[15]/HP(Ωsd2H,1000)[1])^(1/(15))-1)
100*((ΩsdH[15]/ΩsdH[1])^(1/(15))-1)
# 1981-1992
time_annual[34:45]
100*((HP(Ωsd2H,1000)[30]/HP(Ωsd2H,1000)[19])^(1/(12))-1)
100*((ΩsdH[30]/ΩsdH[19])^(1/(12))-1)
# 1993-2008
time_annual[35:end]
100*((HP(Ωsd2H,1000)[end]/HP(Ωsd2H,1000)[31])^(1/(27))-1)
100*((ΩsdH[end]/ΩsdH[31])^(1/(27))-1)

## markup
# 1963-1977
time_annual[16:30]
100*((kysd2H[15]/kysd2H[1])^(1/(15))-1)
100*((kysdH[15]/kysdH[1])^(1/(15))-1)
# 1981-1992
time_annual[34:45]
100*((kysd2H[30]/kysd2H[19])^(1/(12))-1)
100*((kysdH[30]/kysdH[19])^(1/(12))-1)
# 1993-2008
time_annual[35:end]
100*((kysd2H[end]/kysd2H[31])^(1/(27))-1)
100*((kysdH[end]/kysdH[31])^(1/(27))-1)

## unemployment
# 1963-1977
time_annual[16:30]
100*((HP(Usd2H,1000)[15]/HP(Usd2H,1000)[1])^(1/(15))-1)
100*((UsdH[15]/UsdH[1])^(1/(15))-1)
# 1981-1992
time_annual[34:45]
100*((HP(Usd2H,1000)[30]/HP(Usd2H,1000)[19])^(1/(12))-1)
100*((UsdH[30]/UsdH[19])^(1/(12))-1)
# 1993-2008
time_annual[35:end]
100*((HP(Usd2H,1000)[end]/HP(Usd2H,1000)[31])^(1/(27))-1)
100*((UsdH[end]/UsdH[31])^(1/(27))-1)

		############### Code Related to Section 5

		data_premium=data_calibration["wage_premium_acemogluautor!A1:D231"]
		clphsg_all=convert(Array{Float64,1},data_premium[2:end-10,2])
		clghsg_all3=convert(Array{Float64,1},data_premium[2:end-10,4])
		
		year_premium=convert(Array{Float64,1},data_premium[2:end-10,1])
		
function pred_fun(variable, dep)
	n=length(variable)
	X=[ones(n) variable]
	β = inv(X'X)*(X'dep)
	return X*β
end

college_prem_pred = [ones(44) μsd2H[1:end-2]]*(inv([ones(44) μsd2H[1:end-2]]'*[ones(44) μsd2H[1:end-2]])*([ones(44)  μsd2H[1:end-2]]'*unique(clphsg_all[1:4:end])))


with_theme(theme_latexfonts()) do       
    f = Figure(size = (500, 300))
    ax = Axis(f[1, 1],
        xlabel = "Period ",xgridstyle = :dash, ygridstyle = :dash,
        title = "College Premium",
    )
    lines!(ax,time_annual[16:end-2],HP(college_prem_pred,1000), label = "Prediction")
    lines!(ax, time_annual[16:end-2],unique(clphsg_all[1:4:end]), label = "College premium")
    axislegend(; position = :lt, nbanks = 2, framecolor = (:grey, 0.95))
    f
    end		

# Inequality data
top10IncUSA=convert(Array{Float64,1},data_annual[12:74,18]) # 148-2010
top1IncUSA=convert(Array{Float64,1},data_annual[12:74,19]) # 1948-2010


	with_theme(theme_latexfonts()) do       
		f = Figure(size = (500, 300))
		ax = Axis(f[1, 1],
			xlabel = "Period ",xgridstyle = :dash, ygridstyle = :dash,
			title = "Top 1 percent",
		)
		lines!(ax,time_annual[16:end],HP(pred_fun(μsd2H[1:end],top1IncUSA[16:end-2]),1000), label = "Prediction")
		lines!(ax, time_annual[16:end],top1IncUSA[16:end-2], label = "Data")
		axislegend(; position = :lt, nbanks = 2, framecolor = (:grey, 0.95))
		f
		end		



		with_theme(theme_latexfonts()) do       
			f = Figure(size = (500, 300))
			ax = Axis(f[1, 1],
				xlabel = "Period ",xgridstyle = :dash, ygridstyle = :dash,
				title = "Top 10 percent",
			)
			lines!(ax,time_annual[16:end],HP(pred_fun(μsd2H[1:end],top10IncUSA[16:end-2]),1000), label = "Prediction")
			lines!(ax, time_annual[16:end],top10IncUSA[16:end-2], label = "Data")
			axislegend(; position = :lt, nbanks = 2, framecolor = (:grey, 0.95))
			f
			end		



		function r_squared(y_actual::Vector, y_predicted::Vector)
			# Ensure the vectors are the same length
			@assert length(y_actual) == length(y_predicted) 
			ss_res = sum((y_actual .- y_predicted).^2)
			ss_tot = sum((y_actual .- mean(y_actual)).^2)
			return 1 - ss_res / ss_tot
		end
		

# CAGR
		100*((pred_fun(Mech[1:end-7],top10IncUSA[16:end-2])[31]/pred_fun(Mech[1:end-7],top10IncUSA[16:end-2])[19])^(1/13) -1)
		100*((pred_fun(Mech[1:end-7],top1IncUSA[16:end-2])[31]/pred_fun(Mech[1:end-7],top1IncUSA[16:end-2])[19])^(1/13) -1)
		100*((pred_fun(Mech[1:end-9],unique(clphsg_all[1:4:end]))[31]/pred_fun(Mech[1:end-9],unique(clphsg_all[1:4:end]))[19])^(1/13) -1)
		100*(((unique(clphsg_all[1:4:end]))[31]/(unique(clphsg_all[1:4:end]))[19])^(1/13) -1)
		
		

#using DataFrames,CSV,JLD2,DelimitedFiles
#dt_steady_state_dechez = DataFrame(years_dt = time_annual[16:end],ls_beablsd = w_share[16:end-10], ls_full=HP(Ωsd2H,1000), ls_auto=ΩsdH, markup_data=μ_data[16:end-9],markup_full=HP(μsd2H,1000),
#markup_auto=μsdH, KY_data= KY[16:end-9], KY_full= kysd2H./12, KY_auto=kysdH./12, Un=un_c[16:end], Un_full= NRU[16:end], Un_auto=UsdH.*100, un_AI=100*HP(Usd2H,1000))
#		CSV.write("/Users/juanjacobo/Documents/Documents/metroeconomica_2025/steady_dechez.csv", string.(dt_steady_state_dechez))


#dt_calibration = DataFrame(years_dt = time_annual[16:end],auto_data = merged_ctry.share_auto95_US[1:end-7], auto_pred=M̂1, Tw_perfect = tsw5H, Tw_pred=exp.(tsw_indiriect))
#		CSV.write("/Users/juanjacobo/Documents/Documents/metroeconomica_2025/dt_calibration.csv", string.(dt_calibration))


#dtpikettydata = DataFrame(years_dt = time_annual[16:end], top10=top10IncUSA[16:end-2], top1=top1IncUSA[16:end-2],
#top10_mu=HP(pred_fun(μsd2H[1:end],top10IncUSA[16:end-2]),1000),
#top1_mu=HP(pred_fun(μsd2H[1:end],top1IncUSA[16:end-2]),1000), top10_auto=pred_fun(Mech[1:end-7],top10IncUSA[16:end-2]),
#top1_auto=pred_fun(Mech[1:end-7],top1IncUSA[16:end-2]))
#		CSV.write("/Users/juanjacobo/Documents/Documents/metroeconomica_2025/dtpikettydata.csv", string.(dtpikettydata))


#dt_colpre = DataFrame(years_dt = time_annual[16:end-2], colprem=unique(clphsg_all[1:4:end]), cp_mu = HP(pred_fun(μsd2H[1:end-2], unique(clphsg_all[1:4:end])),1000),
#dp_auto = pred_fun(Mech[1:end-9],unique(clphsg_all[1:4:end])))
#				CSV.write("/Users/juanjacobo/Documents/Documents/metroeconomica_2025/dt_colpre.csv", string.(dt_colpre))
		