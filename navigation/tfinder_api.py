# Copyright (c) 2023 Minniti Julien

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of TFinder and associated documentation files, to deal
# in TFinder without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of TFinder, and to permit persons to whom TFinder is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of TFinder.

# TFINDER IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH TFINDER OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import streamlit as st


def tfinder_api():
    st.divider()
    st.markdown("<h3 style='text-align: center; color: black;'>Documentation for TFinder python package</h1>",
                unsafe_allow_html=True)
    st.markdown('How to use:')

    st.divider()
    st.markdown("<h3 style='text-align: center; color: black;'>DNA Region Extractor</h1>",
                unsafe_allow_html=True)
    with st.expander('Analyse gene availability'):
        st.text('Some genes do not have the same name in different species. It can also happen that the gene ID is incorrect.')
        analyse_gene = 'NCBIdna(gene_id).analyse_gene()'
        st.code(analyse_gene)
        r"""Insert containers separated into tabs.

                Inserts a number of multi-element containers as tabs.
                Tabs are a navigational element that allows users to easily
                move between groups of related content.

                To add elements to the returned containers, you can use "with" notation
                (preferred) or just call methods directly on the returned object. See
                examples below.

                .. warning::
                    All the content of every tab is always sent to and rendered on the frontend.
                    Conditional rendering is currently not supported.

                Parameters
                ----------
                tabs : list of strings
                    Creates a tab for each string in the list. The first tab is selected by default.
                    The string is used as the name of the tab and can optionally contain Markdown,
                    supporting the following elements: Bold, Italics, Strikethroughs, Inline Code,
                    Emojis, and Links.

                    This also supports:

                    * Emoji shortcodes, such as ``:+1:``  and ``:sunglasses:``.
                      For a list of all supported codes,
                      see https://share.streamlit.io/streamlit/emoji-shortcodes.

                    * LaTeX expressions, by wrapping them in "$" or "$$" (the "$$"
                      must be on their own lines). Supported LaTeX functions are listed
                      at https://katex.org/docs/supported.html.

                    * Colored text, using the syntax ``:color[text to be colored]``,
                      where ``color`` needs to be replaced with any of the following
                      supported colors: blue, green, orange, red, violet, gray/grey, rainbow.

                    Unsupported elements are unwrapped so only their children (text contents) render.
                    Display unsupported elements as literal characters by
                    backslash-escaping them. E.g. ``1\. Not an ordered list``.

                Returns
                -------
                list of containers
                    A list of container objects.

                Examples
                --------
                You can use `with` notation to insert any element into a tab:

                >>> import streamlit as st
                >>>
                >>> tab1, tab2, tab3 = st.tabs(["Cat", "Dog", "Owl"])
                >>>
                >>> with tab1:
                ...    st.header("A cat")
                ...    st.image("https://static.streamlit.io/examples/cat.jpg", width=200)
                ...
                >>> with tab2:
                ...    st.header("A dog")
                ...    st.image("https://static.streamlit.io/examples/dog.jpg", width=200)
                ...
                >>> with tab3:
                ...    st.header("An owl")
                ...    st.image("https://static.streamlit.io/examples/owl.jpg", width=200)

                .. output ::
                    https://doc-tabs1.streamlit.app/
                    height: 620px

                Or you can just call methods directly in the returned objects:

                >>> import streamlit as st
                >>> import numpy as np
                >>>
                >>> tab1, tab2 = st.tabs(["📈 Chart", "🗃 Data"])
                >>> data = np.random.randn(10, 1)
                >>>
                >>> tab1.subheader("A tab with a chart")
                >>> tab1.line_chart(data)
                >>>
                >>> tab2.subheader("A tab with the data")
                >>> tab2.write(data)


                .. output ::
                    https://doc-tabs2.streamlit.app/
                    height: 700px

                """

