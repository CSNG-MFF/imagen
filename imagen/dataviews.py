import param
import bisect

from imagen.sheetcoords import SheetCoordinateSystem

try:     from collections import OrderedDict
except:  OrderedDict = None


class DataView(param.Parameterized):
    """
    A dataview associates data with the appropriate bounds.
    """

    timestamp = param.Number(default=None, doc=
        """ The initial timestamp. May be used to check for staleness.""")

    bounds = param.Parameter(default=None, doc=
        """ The bounds of the two dimensional coordinate system in which the data resides.""")

    cyclic_interval = param.Number(default=None, allow_None=True, doc=
        """ If None, the data is in a non-cyclic dimension.  Otherwise, the data
            is in a cyclic dimension over the given interval.""")

    roi = param.Parameter(default=None, doc=
        """ The bounds of the region of interest (if any). Any attempt to access
            data outside this region will generate a warning. """)

    labels = param.Dict(default={}, doc=
       """ Additional labels to be associated with the Dataview. """)

    style = param.Dict(default={}, doc=
       """ Metadata hints for visualization purposes.""")

    def __init__(self, bounds, **kwargs):
        super(DataView,self).__init__(bounds=bounds, **kwargs)
        self._data = None

    def record(self, data):
        """
        Records the data into the DataView.
        """
        self._data = data

    def sample(self, coords):
        """
        Return the data at a specified coordinate point over the
        appropriate time slice.
        """
        return data[self._coords2matrixidx(coords, data)]

    def _coords2matrixidx(self, coords):
        raise NotImplementedError

    def view(self):
        """
        Return the requested view as a (data, bbox) tuple.  Provided
        for backward compatibility with the original Topographica
        SheetView model. It is now easier to access the data and
        bounds attributes directly.
        """
        return (self._data, self.bounds)


class IndexedView(DataView):

    indexed_feature = param.String(default='time', doc="""
          The name of the feature used to index the DataView. By
          default this feature is the temporal dimension.""")

    indexed_ascending = param.Boolean(default=True, doc="""
          Whether the indexing feature is to be sorted when returned
          from a slice, regardless of the order of recording. This
          makes sense for features such as time which always has a
          clear ordering.""")

    def __init__(self, bounds, **kwargs):
        super(IndexedView,self).__init__(bounds=bounds, **kwargs)
        self.map_type = OrderedDict if OrderedDict else list
        self._data = []               # Sorted stack if indexed_ascending
        self._index_labels = []       # All the index_labels

    def __getitem__(self, indexslice):
        """
        DataViews support slicing the available data over a time
        interval.  Supports the usual slicing semantics with an
        inclusive lower bound and an exclusive upper bound.
        """
        if not isinstance(indexslice, slice):
            matched_index = (self._bisect_index(indexslice) if self.indexed_ascending else
                             self._index_labels.index(indexslice))
            return self._data[matched_index]

        start, stop, step = indexslice.start, indexslice.stop, indexslice.step

        if self.indexed_ascending:
            return  self._bisect_slice(start, stop, step)
        else:
            (min_val, max_val) = (start, stop) if (start < stop) else (stop, start)
            return self.map_type([pair for pair in zip(self._index_labels, self._data)
                                  if (min_val < l < max_val)])

    def record(self, data, index_value):
        """
        Records the data with the given index_value into the DataView.
        """
        if self.indexed_ascending:
            self._bisect_insert(data, index_value)
        else:
            self._data.append(data.copy())
            self._index_labels.append(index_value)

    def sample(self, coords):
        """
        Return the data at a specified coordinate point over the
        appropriate time slice.
        """
        return [data[self._coords2matrixidx(coords, data)] for data in self._data]

    def _coords2matrixidx(self, coords):
        raise NotImplementedError

    def view(self):
        """
        Return the requested view as a (data, bbox) tuple.  Provided
        for backward compatibility with the original Topographica
        SheetView model. It is now easier to access the data and
        bounds attributes directly.
        """
        return (self._data[-1], self.bounds)

    def __len__(self):
        return len(self._index_labels) if hasattr(self, '_index_labels') else 1

    #========================================================#
    # Bisect module optimizes operations over a sorted index #
    #========================================================#

    def _bisect_index(self, index_value):
        """
        Locate the leftmost value of the stack matching the index_value.
        If not found, raises a ValueError. See Python documentation on
        bisect module for more information about this function.
        """
        i = bisect.bisect_left(self._index_labels, index_value)
        if i != len(self._index_labels) and self._index_labels[i] == index_value: return i
        raise ValueError

    def _bisect_slice(self, start, stop, step):

        start_ind = None if (start is None) else bisect.bisect_left(self._index_labels, start)
        stop_ind = None if (stop is None) else bisect.bisect_left(self._index_labels, stop)
        indexslice = slice(start_ind, stop_ind, step)
        return self.map_type(zip(self._index_labels[indexslice],
                                 self._data[indexslice]))

    def _bisect_insert(self, data, index_value):
        insert_index = bisect.bisect_left(self._index_labels, index_value)
        self._data.insert(insert_index, data.copy())
        self._index_labels.insert(insert_index, index_value)


class Cartesian2D(DataView):
    """
    A dataview for data situated in a 2-dimensional Cartesian
    coordinate system.
    """

    def _coords2matrixidx(self, coords, arr):
        (l,b,r,t) = self.bounds.lbrt()
        (dim1, dim2) = arr.shape
        xdensity = dim1 / (r-l)
        ydensity = dim2 / (t-b)
        return SheetCoordinateSystem(self.bounds, xdensity, ydensity).sheet2matrixidx(*coords)


# Suggested convention: Indexable DataViews can have an x suffix?

# Could make IndexedView into a mix-in class instead...

class Cartesian2Dx(IndexedView):
    _coords2matrixidx = Cartesian2D._coords2matrixidx
