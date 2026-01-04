## ABOUT
See <a href="https://gotankersley.github.io/sloi/2025/11/01/entropic-transforms.html">this post</a> for more information.

<b>Note:</b> This method is only an approximation of entropic order, and not exact.

## PERFORMANCE
This is an optimized version which uses the excellent and amazing <a href="https://flintlib.org/" target="_blank">FLINT Library</a> for calculating large numbers.

There are two main transforms:
1. Near Entropic Order - Implements the ranking using a Restricted Growth Function (RGF)<br>
This allows for sequences in the range of about [100,000 - 200,000], using an alphabet of 16 symbols. (i.e. Hexadecimal).
This is around 50 - 100 kb of input message.

2. Nearer Entropic Order - Implements ranking using a novel method for set partions.<br>
This allows for sequences in the range of about [1,000 - 2,000], using an alphabet of 16 symbols.


## GAINS
Unfortunately, while the transformed message has less entropy on average than the input message, it is often very close, with very small gains. 


## REFERENCES
See <a target="_blank" href="https://ieeexplore.ieee.org/document/9360545">Three Representations for Set Partitions</a> for a mathematical framework of how this works.
Also see <a target="_blank" href="https://www.researchgate.net/publication/377851895_Lexicographic_unranking_algorithms_for_the_Twelvefold_Way">Lexicographic_unranking_algorithms_for_the_Twelvefold_Way</a> as well.