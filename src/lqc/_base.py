"""Internal base for the four labelled accumulator classes."""


class _LabelledStat:
    def __init__(self, label=''):
        if not isinstance(label, str):
            raise TypeError('label should be string.')
        self.label = label

    def __radd__(self, other):
        if other == 0:
            return self
        return self.__add__(other)

    def _require_same_type(self, other):
        if not isinstance(other, type(self)):
            raise TypeError('wrong object to add')
        return other