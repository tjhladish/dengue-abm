# dengue
---
This stochastic, discrete-time, agent-based model of dengue transmission explicitly represents a synthetic population of humans and mosquitoes with which transmission of all four Dengue serotypes can be simulated.

## Daily Transmission Cycle
This ABM is structured around a daily cycle of human-to-mosquito and mosquito-to-human transmission.
Uninfected mosqutioes can acquire dengue after biting infectious humans and can further spread the virus to susceptible humans.
Mosquito infection follows an SEI model of infection progressing between **S**usceptible, **E**xposed, and **I**nfectious.
Human infection similarly follows an SEIR model of infection progressing between **S**usceptible, **E**xposed, **I**nfectious, and **R**esistant. Moreover, human infections can develop varying levels of severity: asymptomatic (**I<sub>A</sub>**), symptomatic (**I<sub>S</sub>**), and withdrawn (**I<sub>W</sub>**).
The model also uses fitted seasonal rain and temperature patterns to reliably determine the synthetic mosquito population size and the extrinsic incubation period (EIP) for mosquito dengue infections.

## Human Synthetic Population
The human synthetic population for this ABM is assumed to have a constant size and age structure in order to reduce computational complexity. Individuals have fixed age, gender, and household assignment.
Within the ABM's daily simulation cycle, humans move between their designated home and a daytime location (either a workplace or a school).
Instead of modeling an aging population, immune histories are annually shifted from younger to older people in order to capture the accumulation of immunity over time which is vital for modeling long-term transmission dynamics.

## Mosquito Synthetic Population
The overall mosquito population size is calculated from location-specific sizes that vary based on seasonal rainfall patterns. Upon infection, a mosquito will be assigned an current age, EIP, and age of death sanmpled from respective distributions.
Mosquito movement can occur between nearby locations weighted by the inverse square distance to a specific location from the origin.
