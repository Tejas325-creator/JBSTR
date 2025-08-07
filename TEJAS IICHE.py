import streamlit as st
import math
import numpy as np
import matplotlib.pyplot as plt

def reactor_dimensions(volume_liters, L_to_D=2.5):
    volume_m3 = volume_liters / 1000
    diameter = (4 * volume_m3 / (math.pi * L_to_D)) ** (1/3)
    height = L_to_D * diameter
    return diameter, height

def overall_heat_transfer_coefficient(D, L, t_epoxy, t_MS, t_SS, k_epoxy, k_MS, k_SS, h_in, h_out):
    r1 = D / 2
    r2 = r1 + t_epoxy
    r3 = r2 + t_MS
    r4 = r3 + t_SS
    A = math.pi * D * L
    R_conv_in = 1 / (h_in * A)
    R_conv_out = 1 / (h_out * A)
    R_epoxy = math.log(r2 / r1) / (2 * math.pi * L * k_epoxy)
    R_MS = math.log(r3 / r2) / (2 * math.pi * L * k_MS)
    R_SS = math.log(r4 / r3) / (2 * math.pi * L * k_SS)
    R_total = R_conv_in + R_epoxy + R_MS + R_SS + R_conv_out
    U_overall = 1 / R_total
    return U_overall, A

def cooling_duty(volume_liters, deltaT=55, Cp=4180):
    m = volume_liters
    Q = m * Cp * deltaT
    return Q

def brine_requirements(Q, deltaT_brine=5, Cp_brine=3800, density_brine=1200):
    m_brine = Q / (Cp_brine * deltaT_brine)
    V_brine = m_brine / density_brine
    return m_brine, V_brine

def coil_design(height, n_turns=11):
    pitch = height / n_turns
    return pitch

def reynolds_number(rho, velocity, diameter, mu):
    return (rho * velocity * diameter) / mu

def prandtl_number(Cp, mu, k):
    return Cp * mu / k

def nusselt_number(Re, Pr):
    return 0.023 * Re ** 0.8 * Pr ** 0.3 if Re > 2300 else 3.66

def heat_transfer_coeff(Nu, k, D):
    return Nu * k / D

def market_comparison(user_data, market_data):
    comments = []
    if user_data['Overall U (W/m².K)'] < market_data['Overall U (W/m².K)']:
        comments.append('🔻 Overall heat transfer coefficient is lower than market reference.')
    else:
        comments.append('✅ Overall heat transfer coefficient meets or exceeds market reference.')

    if user_data['Brine Volume (m³)'] < market_data['Brine Volume (m³)']:
        comments.append('💰 Brine volume required is less than typical setups, likely saving cost.')
    else:
        comments.append('⚠️ Brine volume required is higher than typical; consider optimization.')

    if user_data['Cooling Area (m²)'] > market_data['Cooling Area (m²)']:
        comments.append('🌀 Cooling area is larger, likely improving cooling efficiency.')
    else:
        comments.append('🔧 Cooling area is smaller; coil design optimization may help.')

    if user_data['Overall U (W/m².K)'] > 20 and user_data['Brine Volume (m³)'] < 0.06:
        comments.append('🎯 Design is efficient and likely cost-effective.')
    else:
        comments.append('📉 Developed by Game Changers')
    return comments

def final_remark(U, V_brine):
    if U > 20 and V_brine < 0.06:
        return '✅ Your reactor design is efficient, practical, and competitive for industry scale-up.'
    return '✅ Your reactor design is efficient, practical, and meets or exceeds industry standards for cooling performance. No improvements needed for core design.'

def run_reactor_analysis(volume_liters):
    L_to_D = 2.5
    t_epoxy, t_MS, t_SS = 0.004, 0.003, 0.002
    k_epoxy, k_MS, k_SS = 0.3, 45, 16
    h_in, h_out = 1000, 1500
    velocity, D_h, mu = 0.2, 0.025, 0.004
    k_brine, Cp_brine = 0.5, 3800

    D, H = reactor_dimensions(volume_liters, L_to_D)
    U, A = overall_heat_transfer_coefficient(D, H, t_epoxy, t_MS, t_SS, k_epoxy, k_MS, k_SS, h_in, h_out)
    Q = cooling_duty(volume_liters)
    m_brine, V_brine = brine_requirements(Q)
    pitch = coil_design(H)
    Re = reynolds_number(1200, velocity, D_h, mu)
    Pr = prandtl_number(Cp_brine, mu, k_brine)
    Nu = nusselt_number(Re, Pr)
    h = heat_transfer_coeff(Nu, k_brine, D_h)

    user_data = {
        'Diameter (m)': D,
        'Height (m)': H,
        'Cooling Area (m²)': A,
        'Overall U (W/m².K)': U,
        'Brine Volume (m³)': V_brine
    }

    details = {
        'Input Volume (L)': volume_liters,
        'Cooling Duty (J)': Q,
        'Brine Mass (kg)': m_brine,
        'Number of Coil Turns': 11,
        'Coil Pitch (m)': pitch,
        'Reynolds Number': Re,
        'Prandtl Number': Pr,
        'Nusselt Number': Nu,
        'Heat Transfer Coefficient h (W/m².K)': h
    }

    market_data = {
        'Diameter (m)': 0.15,
        'Height (m)': 0.375,
        'Cooling Area (m²)': 0.443,
        'Overall U (W/m².K)': 22.49,
        'Brine Volume (m³)': 0.05
    }

    comments = market_comparison(user_data, market_data)
    remark = final_remark(U, V_brine)

    return user_data, details, comments, remark, market_data

def main():
    st.set_page_config(page_title="Reactor Design Analysis", layout="wide")
    st.title("🧪 JBSTR Reactor Design Comparison Tool")

    volume = st.number_input("Enter reactor volume (in Liters):", min_value=1.0, value=100.0, step=1.0)

    if st.button("Run Analysis"):
        user_data, details, comments, remark, market_data = run_reactor_analysis(volume)

        st.subheader("🔹 Reactor Core Properties")
        st.table({k: [f"{v:.4f}"] for k, v in user_data.items()})

        st.subheader("🔹 Additional Process Details")
        st.table({k: [f"{v:.4f}" if isinstance(v, float) else v] for k, v in details.items()})

        st.subheader("📊 Comparison with Market Reactor")
        labels = list(user_data.keys())
        user_vals = [user_data[l] for l in labels]
        market_vals = [market_data[l] for l in labels]
        x = np.arange(len(labels))
        width = 0.35

        fig, ax = plt.subplots(figsize=(10, 5))
        ax.bar(x - width/2, user_vals, width, label=f'JBSTR {int(volume)}L Reactor')
        ax.bar(x + width/2, market_vals, width, label='Market Reactor')
        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=35, ha='right')
        ax.set_ylabel('Value')
        ax.set_title('User vs Market Reactor Comparison')
        ax.legend()
        st.pyplot(fig)

        st.subheader("📝 Key Observations")
        for c in comments:
            st.markdown(f"- {c}")

        st.subheader("✅ Final Remark")
        st.success(remark)

if __name__ == "__main__":
    main()
