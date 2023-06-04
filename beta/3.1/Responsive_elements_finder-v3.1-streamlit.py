import streamlit as st

#StreamLit

st.title('Responsive Elements Finder')

# Section "Promoter finder"
st.sidebar.markdown("# Promoter Finder")

# Gene ID entry
gene_id = st.sidebar.text_area("Gene ID:")

# Species selection
species = st.sidebar.selectbox("Species:", ["Human", "Mouse", "Rat"])

# Upstream/downstream entry
upstream = st.sidebar.number_input("Upstream (bp):", value=2000)
downstream = st.sidebar.number_input("Downstream (bp):", value=500)

# Search
if st.sidebar.button("Find promoter"):
    # Perform promoter search here
    pass

# Section "Responsive Elements finder"
st.sidebar.markdown("# Responsive Elements Finder")

# RE entry
responsive_element = st.sidebar.text_input("Responsive element (IUPAC authorize):", value="ATGCCGTA")

# TIS entry
tis = st.sidebar.number_input("Transcription Initiation Site (bp)", value=0)
st.sidebar.markdown("(distance from start of promoter or 'Upstream' if you use the Promoter Finder)")

# Threshold
threshold = st.sidebar.number_input("Threshold (%)", value=80)

# Find RE
if st.sidebar.button("Find responsive elements"):
    # Perform responsive elements search here
    pass

# Section "Status"
st.sidebar.markdown("# Status")

# Configure grid weights
st.beta_set_page_config(layout="wide")

# Section "Promoter finder"
st.subheader("Promoter Finder")

# Gene ID output
st.text_area("Gene ID:", value=gene_id, height=5)

# Species output
st.text("Species: {}".format(species))

# Upstream/downstream output
col1, col2 = st.beta_columns(2)
col1.text("Upstream (bp): {}".format(upstream))
col2.text("Downstream (bp): {}".format(downstream))

# Section "Responsive Elements finder"
st.subheader("Responsive Elements Finder")

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
