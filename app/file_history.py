"""One undo / redo history for a file of several records (chromosomes).

Each record keeps its own undo steps in its Document (Rust). This log only remembers in which order they were made, so that
Edit > Undo undoes the last edit of the whole file, whichever chromosome it was in (the window then shows that chromosome), and
the History panel lists the whole file's changes oldest first.

It learns from what the editors report after each change: (record, undo steps, redo steps). Nothing here needs Qt.
"""


class FileHistory:
    def __init__(self):
        self.undo_log = []     # record of each undo step of the file, oldest first
        self.redo_log = []     # record of each undone step; the last one is the next to redo
        self.seen = {}         # record -> (undo steps, redo steps) at the last report

    def reset(self):
        self.undo_log, self.redo_log, self.seen = [], [], {}

    def forget(self, record):
        """The record was removed from the file: its steps leave the file's history."""
        self.undo_log = [r for r in self.undo_log if r != record]
        self.redo_log = [r for r in self.redo_log if r != record]
        self.seen.pop(record, None)

    def observe(self, record, undo_steps, redo_steps):
        """The record now has this many undo / redo steps. The first report of a record is only a starting point."""
        before = self.seen.get(record)
        self.seen[record] = (undo_steps, redo_steps)
        if before is None:
            return
        grew = undo_steps - before[0]
        if grew > 0 and redo_steps == before[1] - grew:     # steps redone
            self._move(self.redo_log, self.undo_log, record, grew)
        elif grew > 0:                                      # new edits: what was undone in that record cannot be redone
            self.undo_log += [record] * grew
            self.redo_log = [r for r in self.redo_log if r != record]
        elif grew < 0:                                      # steps undone
            self._move(self.undo_log, self.redo_log, record, -grew)

    @staticmethod
    def _move(source, target, record, count):
        """Move the last `count` entries of `record` from one list to the other (newest last)."""
        for _ in range(count):
            for i in range(len(source) - 1, -1, -1):
                if source[i] == record:
                    target.append(source.pop(i))
                    break

    def next_undo(self):
        return self.undo_log[-1] if self.undo_log else None

    def next_redo(self):
        return self.redo_log[-1] if self.redo_log else None

    def steps(self, histories):
        """The whole file's steps as ([(record, step)] oldest first, [(record, step)] in redo order), taking each record's own
        steps from `histories[record] = (undo list, redo list)` as Document.history() gives them."""
        def merge(log, lists):
            taken, out = {}, []
            for record in log:
                i = taken.get(record, 0)
                taken[record] = i + 1
                steps = lists(record)
                if i < len(steps):
                    out.append((record, steps[i]))
            return out
        undo = merge(self.undo_log, lambda r: histories.get(r, ([], []))[0])
        redo = merge(self.redo_log[::-1], lambda r: histories.get(r, ([], []))[1])
        return undo, redo
