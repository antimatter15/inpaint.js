inpaint.js
==========

This is mostly a port of the Telea Inpainting algorithm from [scikit-image](https://github.com/chintak/scikit-image/blob/inpaint/skimage/filter/inpaint.py). It uses the fast marching method to fill in particular regions of an image based on the surrounding area. The particular details of the implementation are described really well in a [blog post by Chintak Sheth](http://chintaksheth.com/2013/07/fast-marching-method/), the person who originally implemented these inpainting routines in python.

It is in some way similar, and in other ways different from the Content Aware Fill feature in Adobe Photoshop, which in case you haven't seen, has some pretty cool demos that will inevitably allow the Stalinist-Skynet to remove dissidents (i.e. all humans) from photographs. 

Content Aware Fill has methods for texture synthesis (I should totally linkify everything here to Wikipedia at some point) so it's more like the Exemplar based approach used by Criminisi. I think I heard a video somewhere that they used the method by Wexler, et al. but I'm not particularly keen on nailing down the exact algorithm.

It's an algorithm that operates at a pretty reasonable speed, and it's also pretty concise, which is nice. It has a dependency on jsfeat, where the extent of that dependency is essentially to have a nice way to wrap a typed array into something that has a tenuous semblance to a multidimensional array.

Realizing that, there's actually a pretty good chance that I'll get rid of that dependency because it'd literally entail nothing more than the replacement of like four lines of code. 

