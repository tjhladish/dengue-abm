death_p = c(
  0.0018, 0.002069514,0.002378646,0.002732982,0.003138823,0.003603251,0.004134191,
  0.004740476,0.005431895,0.006219231,0.007114273,0.008129798,0.009279508,
  0.01057792,0.01204018,0.01368180,0.01551828,0.01756467,0.01983494,0.02234135,
  0.02509359,0.02809799,0.03135655,0.03486617,0.03861782,0.04259607,0.04677877,
  0.05113718,0.05563645,0.06023658,0.0648937,0.0695617,0.07419404,0.07874562,
  0.08317441,0.08744304,0.09151983,0.09537952,0.09900354,0.1023799,0.1055030,
  0.1083723,0.1109923,0.1133714,0.1155206,0.1174531,0.1191837,0.1207277,0.1221007,
  0.1233179,0.1243943,0.1253438,0.1261799,0.1269147,0.1275594,0.1281243,0.1286187,
  0.1290510,0.1294285, 1 #0.1297580, 0.1300453
)

survives_p = 1 - death_p
survives_from_birth = cumprod(survives_p)
frac_alive = c(1, head(survives_from_birth,-1)) # for daily mu calc

bite_age_cdf = cumsum(frac_alive/sum(frac_alive))

death_age_pdf = c(death_p[1], head(survives_to_p,-1)*tail(death_p,-1))
death_age_cdf = cumsum(death_age_pdf)

# for biting calculations:
# biting age (without mod) draw against biting_age_cdf
# death age: take biting age, draw from death_age_cdf[biting_age-1] to 1 (0 if biting age = 1)

# for daily death mu:
# M*(1-eff) = M'
# M = sum(frac_alive), M' = sum(frac_alive * (1-mu)^((1:length(frac_alive))-1))

uniroot(function(mu, coeffs, eff){
  sum(coeffs)*(1-eff) - sum(coeffs * (1-mu)^((1:length(coeffs))-1))
}, c(0,1), coeffs=frac_alive, eff=0.8)