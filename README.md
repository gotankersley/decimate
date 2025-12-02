## ABOUT
See <a href="https://gotankersley.github.io/sloi/2025/11/01/entropic-transforms.html">this post</a> for more information.

<b>Note:</b> This method is only an approximation of entropic order, and not exact.

## PERFORMANCE
This is an optimized version which uses the excellent and amazing <a href="https://flintlib.org/" target="_blank">FLINT Library</a> for calculating large numbers.

This allows for sequences in the range of about [100,000 - 200,000], using an alphabet of 16 symbols. (i.e. Hexadecimal).
This is around 50 - 100 kb of input message.

## GAINS
Unfortunately, while the transformed message has less entropy on average than the input message, it is often very close, with very small gains. 


## REFERENCES
See <a href="https://ieeexplore.ieee.org/document/9360545">Three Representations for Set Partitions</a> for a mathematical framework of how this works.