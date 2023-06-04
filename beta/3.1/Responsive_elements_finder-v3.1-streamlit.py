import streamlit as st

#StreamLit

st.title('Responsive Elements Finder')

import streamlit as st

# Configure grid layout
col1, col2 = st.beta_columns(2)

# Section "Promoter finder"
with col1:
    st.markdown("# Promoter Finder")

    # Gene ID entry
    gene_id = st.text_area("Gene ID:")

    # Species selection
    species = st.selectbox("Species:", ["Human", "Mouse", "Rat"])

    # Upstream/downstream entry
    upstream = st.number_input("Upstream (bp):", value=2000)
    downstream = st.number_input("Downstream (bp):", value=500)

    # Search
    if st.button("Find promoter"):
        # Perform promoter search here
        pass

# Section "Responsive Elements finder"
with col2:
    st.markdown("# Responsive Elements Finder")

    # RE entry
    responsive_element = st.text_input("Responsive element (IUPAC authorize):", value="ATGCCGTA")

    # TIS entry
    tis = st.number_input("Transcription Initiation Site (bp)", value=0)
    st.markdown("(distance from start of promoter or 'Upstream' if you use the Promoter Finder)")

    # Threshold
    threshold = st.number_input("Threshold (%)", value=80)

    # Find RE
    if st.button("Find responsive elements"):
        # Perform responsive elements search here
        pass

# Section "Status"
st.markdown("# Status")

# Section "Promoter"
st.subheader("Promoter")

# Gene ID output
st.text_area("Gene ID:", value=gene_id, height=5)

# Species output
st.text("Species: {}".format(species))

# Upstream/downstream output
col1, col2 = st.beta_columns(2)
col1.text("Upstream (bp): {}".format(upstream))
col2.text("Downstream (bp): {}".format(downstream))

# Section "Responsive elements"
st.subheader("Responsive Elements")

# Responsive element output
st.text("Responsive element (IUPAC authorize): {}".format(responsive_element))

# TIS output
st.text("Transcription Initiation Site (bp): {}".format(tis))

# Threshold output
st.text("Threshold (%): {}".format(threshold))

# Section "Promoter"
st.subheader("Promoter")

# Promoter output
st.text_area("Promoter:", value="", height=10)

# Section "Responsive elements"
st.subheader("Responsive elements")

# Responsive elements output
st.text_area("Responsive elements:", value="", height=16)

# Section "Status"
st.subheader("Status")

# Status output
st.text("Status: ")

# Help & Credit
st.sidebar.markdown("# Help & Credit")
st.sidebar.button("How to use")
st.sidebar.button("Github by MINNITI Julien")

