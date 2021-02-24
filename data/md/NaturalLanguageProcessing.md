## CRAN Task View: Natural Language Processing

  ----------------- ------------------------------------------------------------------------------------------------------------------------------------------
  **Maintainer:**   Fridolin Wild, Performance Augmentation Lab (PAL, Department of Computing and Communications Technologies, Oxford Brookes University, UK
  **Contact:**      wild at brookes.ac.uk
  **Version:**      2017-11-29
  **URL:**          <https://CRAN.R-project.org/view=NaturalLanguageProcessing>
  ----------------- ------------------------------------------------------------------------------------------------------------------------------------------

<div>

Natural language processing has come a long way since its foundations
were laid in the 1940s and 50s (for an introduction see, e.g., Jurafsky
and Martin (2008): Speech and Language Processing, Pearson Prentice
Hall). This CRAN task view collects relevant R packages that support
computational linguists in conducting analysis of speech and language on
a variety of levels - setting focus on words, syntax, semantics, and
pragmatics.

In recent years, we have elaborated a framework to be used in packages
dealing with the processing of written material: the package
[tm](../packages/tm/index.html). Extension packages in this area are
highly recommended to interface with tm\'s basic routines and useRs are
cordially invited to join in the discussion on further developments of
this framework package. To get into natural language processing, the
[cRunch service](http://cRunch.kmi.open.ac.uk) and
[tutorials](http://cRunch.kmi.open.ac.uk/w/index.php/Tutorials) may be
helpful.

#### Frameworks:

-   [tm](../packages/tm/index.html) provides a comprehensive text mining
    framework for R. The [Journal of Statistical
    Software](http://www.jstatsoft.org/) article [Text Mining
    Infrastructure in R](http://www.jstatsoft.org/v25/i05/) gives a
    detailed overview and presents techniques for count-based analysis
    methods, text clustering, text classification and string kernels.
-   [tm.plugin.dc](../packages/tm.plugin.dc/index.html) allows for
    distributing corpora across storage devices (local files or Hadoop
    Distributed File System).
-   [tm.plugin.mail](../packages/tm.plugin.mail/index.html) helps with
    importing mail messages from archive files such as used in
    Thunderbird (mbox, eml).
-   [tm.plugin.alceste](../packages/tm.plugin.alceste/index.html) allows
    importing text corpora written in a file in the Alceste format.
-   [tm.plugin.factiva](../packages/tm.plugin.factiva/index.html),
    [tm.plugin.lexisnexis](../packages/tm.plugin.lexisnexis/index.html),
    [tm.plugin.europresse](../packages/tm.plugin.europresse/index.html)
    allow importing press and Web corpora from (respectively) Dow Jones
    Factiva, LexisNexis, and Europresse.
-   [tm.plugin.webmining](../packages/tm.plugin.webmining/index.html)
    allow importing news feeds in XML (RSS, ATOM) and JSON formats.
    Currently, the following feeds are implemented: Google Blog Search,
    Google Finance, Google News, NYTimes Article Search, Reuters News
    Feed, Yahoo Finance, and Yahoo Inplay.
-   [RcmdrPlugin.temis](../packages/RcmdrPlugin.temis/index.html) is an
    Rcommander plug-in providing an integrated solution to perform a
    series of text mining tasks such as importing and cleaning a corpus,
    and analyses like terms and documents counts, vocabulary tables,
    terms co-occurrences and documents similarity measures, time series
    analysis, correspondence analysis and hierarchical clustering.
-   [openNLP](../packages/openNLP/index.html) provides an R interface to
    [OpenNLP](http://opennlp.sourceforge.net/) , a collection of natural
    language processing tools including a sentence detector, tokenizer,
    pos-tagger, shallow and full syntactic parser, and named-entity
    detector, using the Maxent Java package for training and using
    maximum entropy models.
-   Trained models for English and Spanish to be used with
    [openNLP](../packages/openNLP/index.html) are available from
    <http://datacube.wu.ac.at/> as packages openNLPmodels.en and
    openNLPmodels.es, respectively.
-   [RWeka](../packages/RWeka/index.html) is a interface to
    [Weka](http://www.cs.waikato.ac.nz/ml/weka/) which is a collection
    of machine learning algorithms for data mining tasks written in
    Java. Especially useful in the context of natural language
    processing is its functionality for tokenization and stemming.
-   [tidytext](../packages/tidytext/index.html) provides means for text
    mining for word processing and sentiment analysis using dplyr,
    ggplot2, and other tidy tools.
-   [monkeylearn](../packages/monkeylearn/index.html) provides a wrapper
    interface to machine learning services on Monkeylearn for text
    analysis, i.e., classification and extraction.
-   [udpipe](../packages/udpipe/index.html) provides
    language-independant tokenization, part of speech tagging,
    lemmatization, dependency parsing, and training of treebank-based
    annotation models.

#### Words (lexical DBs, keyword extraction, string manipulation, stemming)

-   R\'s base package already provides a rich set of character
    manipulation routines. See
    `help.search(keyword = "character", package = "base")` for more
    information on these capabilities.
-   [wordnet](../packages/wordnet/index.html) provides an R interface to
    [WordNet](http://wordnet.princeton.edu/) , a large lexical database
    of English.
-   [RKEA](../packages/RKEA/index.html) provides an R interface to
    [KEA](http://www.nzdl.org/Kea/) (Version 5.0). KEA (for Keyphrase
    Extraction Algorithm) allows for extracting keyphrases from text
    documents. It can be either used for free indexing or for indexing
    with a controlled vocabulary.
-   [gsubfn](../packages/gsubfn/index.html) can be used for certain
    parsing tasks such as extracting words from strings by content
    rather than by delimiters. `demo("gsubfn-gries")` shows an example
    of this in a natural language processing context.
-   [textreuse](../packages/textreuse/index.html) provides a set of
    tools for measuring similarity among documents and helps with
    detecting passages which have been reused. The package implements
    shingled n-gram, skip n-gram, and other tokenizers;
    similarity/dissimilarity functions; pairwise comparisons; minhash
    and locality sensitive hashing algorithms; and a version of the
    Smith-Waterman local alignment algorithm suitable for natural
    language.
-   [boilerpipeR](../packages/boilerpipeR/index.html) helps with the
    extraction and sanitizing of text content from HTML files: removal
    of ads, sidebars, and headers using the boilerpipe Java library.
-   [tau](../packages/tau/index.html) contains basic string manipulation
    and analysis routines needed in text processing such as dealing with
    character encoding, language, pattern counting, and tokenization.
-   [SnowballC](../packages/SnowballC/index.html) provides exactly the
    same API as Rstem, but uses a slightly different design of the C
    libstemmer library from the Snowball project. It also supports two
    more languages.
-   [stringi](../packages/stringi/index.html) provides R language
    wrappers to the International Components for Unicode (ICU) library
    and allows for: conversion of text encodings, string searching and
    collation in any locale, Unicode normalization of text, handling
    texts with mixed reading direction (e.g., left to right and right to
    left), and text boundary analysis (for tokenizing on different
    aggregation levels or to identify suitable line wrapping locations).
-   [stringdist](../packages/stringdist/index.html) implements an
    approximate string matching version of R\'s native \'match\'
    function. It can calculate various string distances based on edits
    (Damerau-Levenshtein, Hamming, Levenshtein, optimal string
    alignment), qgrams (q-gram, cosine, jaccard distance) or heuristic
    metrics (Jaro, Jaro-Winkler). An implementation of soundex is
    provided as well. Distances can be computed between character
    vectors while taking proper care of encoding or between integer
    vectors representing generic sequences.
-   [[Rstem]{.Ohat}](http://www.Omegahat.net/Rstem/) (available from
    Omegahat) is an alternative interface to a C version of Porter\'s
    word stemming algorithm.
-   [KoNLP](../packages/KoNLP/index.html) provides a collection of
    conversion routines (e.g. Hangul to Jamos), stemming, and part of
    speech tagging through interfacing with the Lucene\'s HanNanum
    analyzer. In version 0.0-8.0, the documentation is sparse and still
    needs some help.
-   [koRpus](../packages/koRpus/index.html) is a diverse collection of
    functions for automatic language detection, hyphenation, several
    indices of lexical diversity (e.g., type token ratio, HD-D/vocd-D,
    MTLD) and readability (e.g., Flesch, SMOG, LIX, Dale-Chall). See the
    [web page](http://reaktanz.de/?c=hacking&s=koRpus) for more
    information.
-   [alineR](../packages/alineR/index.html) helps calculate the phonetic
    distance between words (the \'ALINE\' distance). The score is based
    on phonetic featuers represented with the Unicode-compliant
    International Phonetic Alphabet (IPA). Parameterized features
    weights are used to determine the optimal alignment and functions
    are provided to estimate optimum values using a genetic algorithm
    and supervised learning.
-   [ore](../packages/ore/index.html) provides an alternative to R\'s
    built-in functionality for handling regular expressions, based on
    the Onigmo Regular Expression Library. Offers first-class compiled
    regex objects, partial matching and function-based substitutions,
    amongst other features. A benchmark comparing results for ore
    functions with stringi and the R base implementation is available
    [[regex-performance]{.GitHub}](https://github.com/jonclayden/regex-performance/).
-   [languageR](../packages/languageR/index.html) provides data sets and
    functions exemplifying statistical methods, and some facilitatory
    utility functions used in the book by R. H. Baayen: \"Analyzing
    Linguistic Data: a Practical Introduction to Statistics Using R\",
    Cambridge University Press, 2008.
-   [zipfR](../packages/zipfR/index.html) offers some statistical models
    for word frequency distributions. The utilities include functions
    for loading, manipulating and visualizing word frequency data and
    vocabulary growth curves. The package also implements several
    statistical models for the distribution of word frequencies in a
    population. (The name of this library derives from the most famous
    word frequency distribution, Zipf\'s law.)
-   [maxent](../packages/maxent/index.html) is an implementation of
    maxinum entropy minimising memory consumption of very large
    data-sets.
-   [wordcloud](../packages/wordcloud/index.html) provides a
    visualisation similar to the famous wordle ones: it horizontally and
    vertically distributes features in a pleasing visualisation with the
    font size scaled by frequency.
-   [hunspell](../packages/hunspell/index.html) is a stemmer and
    spell-checker library designed for languages with rich morphology
    and complex word compounding or character encoding. The package can
    check and analyze individual words as well as search for incorrect
    words within a text, latex or (R package) manual document.
-   [phonics](../packages/phonics/index.html) provides a collection of
    phonetic algorithms including Soundex, Metaphone, NYSIIS,
    Caverphone, and others.
-   [tesseract](../packages/tesseract/index.html) is an OCR engine with
    unicode (UTF-8) support that can recognize over 100 languages out of
    the box.
-   [mscsweblm4r](../packages/mscsweblm4r/index.html) provides an
    interface to the Microsoft Cognitive Services Web Language Model API
    and can be used to calculate the probability for a sequence of words
    to appear together, the conditional probability that a specific word
    will follow an existing sequence of words, get the list of words
    (completions) most likely to follow a given sequence of words, and
    insert spaces into a string of words adjoined together without any
    spaces (hashtags, URLs, etc.).
-   [mscstexta4r](../packages/mscstexta4r/index.html) provides an
    interface to the Microsoft Cognitive Services Text Analytics API and
    can be used to perform sentiment analysis, topic detection, language
    detection, and key phrase extraction.
-   [tokenizers](../packages/tokenizers/index.html) helps split text
    into tokens, supporting shingled n-grams, skip n-grams, words, word
    stems, sentences, paragraphs, characters, lines, and regular
    expressions.

#### Semantics:

-   [lsa](../packages/lsa/index.html) provides routines for performing a
    latent semantic analysis with R. The basic idea of latent semantic
    analysis (LSA) is, that text do have a higher order (=latent
    semantic) structure which, however, is obscured by word usage (e.g.
    through the use of synonyms or polysemy). By using conceptual
    indices that are derived statistically via a truncated singular
    value decomposition (a two-mode factor analysis) over a given
    document-term matrix, this variability problem can be overcome. The
    article [Investigating Unstructured Texts with Latent Semantic
    Analysis](http://www.springerlink.com/content/g7u377132gq5623g/)
    gives a detailed overview and demonstrates the use of the package
    with examples from the are of technology-enhanced learning.
-   [topicmodels](../packages/topicmodels/index.html) provides an
    interface to the C code for Latent Dirichlet Allocation (LDA) models
    and Correlated Topics Models (CTM) by David M. Blei and co-authors
    and the C++ code for fitting LDA models using Gibbs sampling by
    Xuan-Hieu Phan and co-authors.
-   [lda](../packages/lda/index.html) implements Latent Dirichlet
    Allocation and related models similar to LSA and topicmodels.
-   [stm](../packages/stm/index.html) (Structural Topic Model)
    implements a topic model derivate that can include document-level
    meta-data. The package also includes tools for model selection,
    visualization, and estimation of topic-covariate regressions.
-   [kernlab](../packages/kernlab/index.html) allows to create and
    compute with string kernels, like full string, spectrum, or bounded
    range string kernels. It can directly use the document format used
    by [tm](../packages/tm/index.html) as input.
-   [skmeans](../packages/skmeans/index.html) helps with clustering
    providing several algorithms for spherical k-means partitioning.
-   [movMF](../packages/movMF/index.html) provides another clustering
    alternative (approximations are fitted with von Mises-Fisher
    distributions of the unit length vectors).
-   [RTextTools](../packages/RTextTools/index.html) is a machine
    learning package for automatic text classification. It implements
    the nine different algorithms (svm, slda, boosting, bagging, rf,
    glmnet, tree, nnet, and maxent) and routines supporting the
    evaluation of accuracy.
-   [textir](../packages/textir/index.html) is a suite of tools for text
    and sentiment mining.
-   [textcat](../packages/textcat/index.html) provides support for
    n-gram based text categorization.
-   [textrank](../packages/textrank/index.html) is an extension of the
    PageRank and allows to summarize text by calculating how sentences
    are related to one another.
-   [corpora](../packages/corpora/index.html) offers utility functions
    for the statistical analysis of corpus frequency data.
-   [text2vec](../packages/text2vec/index.html) provides tools for text
    vectorization, topic modeling (LDA, LSA), word embeddings (GloVe),
    and similarities.

#### Pragmatics:

-   [qdap](../packages/qdap/index.html) helps with quantitative
    discourse analysis of transcripts.
-   [quanteda](../packages/quanteda/index.html) supports quantitative
    analysis of textual data.

#### Evaluation:

-   [rel](../packages/rel/index.html) implements a variety of
    performance metrics for assessing prediction quality beyond
    precision and recall, including point estimates with confidence
    intervals for Bennett et al.\'s S, Cohen\'s kappa, Conger\'s kappa,
    Fleiss\' kappa, Gwet\'s AC, intraclass correlation coefficients,
    Krippendorff\'s alpha, Scott\'s pi, the standard error of
    measurement, and weighted kappa.

#### Corpora:

-   [gutenbergr](../packages/gutenbergr/index.html) allows downloading
    and processing public domain works in the Project Gutenberg
    collection. Includes metadata for all Project Gutenberg works, so
    that they can be searched and retrieved.

</div>

### CRAN packages:

-   [alineR](../packages/alineR/index.html)
-   [boilerpipeR](../packages/boilerpipeR/index.html)
-   [corpora](../packages/corpora/index.html)
-   [gsubfn](../packages/gsubfn/index.html)
-   [gutenbergr](../packages/gutenbergr/index.html)
-   [hunspell](../packages/hunspell/index.html)
-   [kernlab](../packages/kernlab/index.html)
-   [KoNLP](../packages/KoNLP/index.html)
-   [koRpus](../packages/koRpus/index.html)
-   [languageR](../packages/languageR/index.html)
-   [lda](../packages/lda/index.html)
-   [lsa](../packages/lsa/index.html)
-   [maxent](../packages/maxent/index.html)
-   [monkeylearn](../packages/monkeylearn/index.html)
-   [movMF](../packages/movMF/index.html)
-   [mscstexta4r](../packages/mscstexta4r/index.html)
-   [mscsweblm4r](../packages/mscsweblm4r/index.html)
-   [openNLP](../packages/openNLP/index.html)
-   [ore](../packages/ore/index.html)
-   [phonics](../packages/phonics/index.html)
-   [phonics](../packages/phonics/index.html)
-   [qdap](../packages/qdap/index.html)
-   [quanteda](../packages/quanteda/index.html)
-   [RcmdrPlugin.temis](../packages/RcmdrPlugin.temis/index.html)
-   [rel](../packages/rel/index.html)
-   [RKEA](../packages/RKEA/index.html)
-   [RTextTools](../packages/RTextTools/index.html)
-   [RWeka](../packages/RWeka/index.html)
-   [skmeans](../packages/skmeans/index.html)
-   [SnowballC](../packages/SnowballC/index.html)
-   [stm](../packages/stm/index.html)
-   [stringdist](../packages/stringdist/index.html)
-   [stringi](../packages/stringi/index.html)
-   [tau](../packages/tau/index.html)
-   [tesseract](../packages/tesseract/index.html)
-   [text2vec](../packages/text2vec/index.html)
-   [textcat](../packages/textcat/index.html)
-   [textir](../packages/textir/index.html)
-   [textrank](../packages/textrank/index.html)
-   [textreuse](../packages/textreuse/index.html)
-   [tidytext](../packages/tidytext/index.html)
-   [tm](../packages/tm/index.html) (core)
-   [tm.plugin.alceste](../packages/tm.plugin.alceste/index.html)
-   [tm.plugin.dc](../packages/tm.plugin.dc/index.html)
-   [tm.plugin.europresse](../packages/tm.plugin.europresse/index.html)
-   [tm.plugin.factiva](../packages/tm.plugin.factiva/index.html)
-   [tm.plugin.lexisnexis](../packages/tm.plugin.lexisnexis/index.html)
-   [tm.plugin.mail](../packages/tm.plugin.mail/index.html)
-   [tm.plugin.webmining](../packages/tm.plugin.webmining/index.html)
-   [tokenizers](../packages/tokenizers/index.html)
-   [topicmodels](../packages/topicmodels/index.html)
-   [udpipe](../packages/udpipe/index.html)
-   [wordcloud](../packages/wordcloud/index.html)
-   [wordnet](../packages/wordnet/index.html)
-   [zipfR](../packages/zipfR/index.html)

### Related links:

-   CRAN Task View: [Cluster](Cluster.html)
-   CRAN Task View: [MachineLearning](MachineLearning.html)
-   Omegahat Package: [[Rstem]{.Ohat}](http://www.Omegahat.net/Rstem/)
-   [The KMi cRunch
    tutorials](http://crunch.kmi.open.ac.uk/w/index.php/Tutorials)
-   [A Gentle Introduction to Statistics for (Computational) Linguists
    (SIGIL)](http://www.stefan-evert.de/SIGIL/)
-   [Stefan Th. Gries (2009): Quantitative Corpus Linguistics with R,
    Routledge.](http://www.routledge.com/books/details/9780415962704/)
-   [ttda: Tools for Textual Data Analysis
    (Deprecated)](http://wwwpeople.unil.ch/jean-pierre.mueller/ttda_-_archives.html)
-   [Corpora and NLP model packages at
    http://datacube.wu.ac.at/](http://datacube.wu.ac.at/)
