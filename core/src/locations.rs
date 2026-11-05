//! How feature locations move when bases are inserted or deleted.
//!
//! Locations use gb-io's `Location` type (the GenBank location grammar as a tree:
//! `complement(join(1..10,20..30))` = Complement(Join([Range, Range]))).
//! Coordinates are 0-based, end-exclusive: GenBank `1..10` is `Range(0, 10)`.
//!
//! Rules (same behaviour as .dna files):
//! - insert strictly inside a range      -> the range grows
//! - insert at or before a range's start -> the range moves right
//! - insert at or after a range's end    -> the range is unchanged
//! - delete overlapping a range          -> the range is trimmed; fully deleted -> removed
//! - a feature whose every part was deleted -> the whole feature is removed (`None`)

use gb_io::seq::{After, Before, Location};

/// New location after inserting `n` bases at position `pos`.
///
/// `match` is Rust's switch statement on steroids: each arm unpacks one variant of the
/// `Location` enum and binds its inner values to names (`s`, `e`, `before`, ...).
pub fn after_insert(loc: &Location, pos: i64, n: i64) -> Location {
    // Small helper closure: a point at or after `pos` moves right by `n`.
    let point = |x: i64| if x >= pos { x + n } else { x };
    match loc {
        Location::Range((s, before), (e, after)) => {
            let (s, e) = (*s, *e); // `*` copies the numbers out of the borrowed tuple
            let (s2, e2) = if pos <= s {
                (s + n, e + n) // insertion before (or at) the start: shift
            } else if pos < e {
                (s, e + n) // insertion strictly inside: grow
            } else {
                (s, e) // insertion after: untouched
            };
            // `.clone()` copies the small Before/After flags (the "<" / ">" fuzzy markers).
            Location::Range((s2, before.clone()), (e2, after.clone()))
        }
        Location::Between(a, b) => Location::Between(point(*a), point(*b)),
        Location::Complement(inner) => Location::Complement(Box::new(after_insert(inner, pos, n))),
        // `.iter().map(...).collect()` builds a new Vec from the transformed children.
        Location::Join(v) => Location::Join(v.iter().map(|l| after_insert(l, pos, n)).collect()),
        Location::Order(v) => Location::Order(v.iter().map(|l| after_insert(l, pos, n)).collect()),
        Location::Bond(v) => Location::Bond(v.iter().map(|l| after_insert(l, pos, n)).collect()),
        Location::OneOf(v) => Location::OneOf(v.iter().map(|l| after_insert(l, pos, n)).collect()),
        // External references point into another record; gaps have no position.
        Location::External(..) | Location::Gap(_) => loc.clone(),
    }
}

/// New location after deleting bases `[start, end)`, or `None` if nothing is left.
///
/// `Option<Location>` is how Rust says "maybe a location": `Some(loc)` or `None`.
pub fn after_delete(loc: &Location, start: i64, end: i64) -> Option<Location> {
    let len = end - start;
    // Where a coordinate lands after the deletion: before it -> same; after it -> shifted
    // left; inside it -> collapses onto `start`.
    let point = |x: i64| {
        if x <= start {
            x
        } else if x >= end {
            x - len
        } else {
            start
        }
    };
    // Shared logic for the list variants: transform each child, drop deleted ones.
    // `filter_map` = map + drop the `None`s in one pass.
    let children = |v: &Vec<Location>| -> Vec<Location> {
        v.iter().filter_map(|l| after_delete(l, start, end)).collect()
    };
    match loc {
        Location::Range((s, before), (e, after)) => {
            let (s2, e2) = (point(*s), point(*e));
            if s2 >= e2 {
                None // the whole range was inside the deletion
            } else {
                Some(Location::Range((s2, before.clone()), (e2, after.clone())))
            }
        }
        Location::Between(a, b) => Some(Location::Between(point(*a), point(*b))),
        // `?` on an Option: if the inner location vanished, the complement vanishes too.
        Location::Complement(inner) => Some(Location::Complement(Box::new(after_delete(inner, start, end)?))),
        Location::Join(v) => collapse(children(v), Location::Join),
        Location::Order(v) => collapse(children(v), Location::Order),
        Location::Bond(v) => non_empty(children(v)).map(Location::Bond),
        Location::OneOf(v) => non_empty(children(v)).map(Location::OneOf),
        Location::External(..) | Location::Gap(_) => Some(loc.clone()),
    }
}

/// A join with no parts left disappears; a join with a single part becomes that part
/// (GenBank `join(5..10)` is legal but ugly).
///
/// `wrap: fn(Vec<Location>) -> Location` = we pass the enum variant itself as a function,
/// e.g. `Location::Join`, so the same code serves Join and Order.
fn collapse(mut parts: Vec<Location>, wrap: fn(Vec<Location>) -> Location) -> Option<Location> {
    match parts.len() {
        0 => None,
        1 => parts.pop(),
        _ => Some(wrap(parts)),
    }
}

fn non_empty(parts: Vec<Location>) -> Option<Vec<Location>> {
    if parts.is_empty() { None } else { Some(parts) }
}

/// The part of `loc` inside a window of `len` bases that starts at `start` on a circle of
/// `n` bases (the window may pass the origin: `start + len > n`), as a location relative
/// to the window start; `None` when nothing of it is inside. This is what Copy across the
/// origin and Rotate use.
///
/// Each part of a join is cut and moved on its own, so the parts stay in their original
/// order (a feature's own sequence is read part after part). A part that straddles the
/// cut splits into two parts, the piece before the cut first. Neighbouring parts that
/// end up touching are merged (`join(1..5,6..10)` reads the same as `1..10`).
pub fn clip_window(loc: &Location, start: i64, len: i64, n: i64) -> Option<Location> {
    let w = Window { start, end: start + len, n };
    let parts = w.parts(loc);
    if parts.is_empty() {
        return None;
    }
    Some(merged(parts, if matches!(loc, Location::Order(_)) { Location::Order } else { Location::Join }))
}

/// The window of `clip_window`: bases `[start, end)` of the circle unrolled (end <= start + n).
struct Window {
    start: i64,
    end: i64,
    n: i64,
}

impl Window {
    /// A range `[s, e)` as 0-2 pieces `(new start, new end, start was cut, end was cut)`:
    /// first the bases before `end - n` (they come after the origin, so they land at the end
    /// of the window), then the bases from `start` up (they land at its start). That is
    /// also the order along the range, so a range keeps reading in the same direction.
    fn pieces(&self, s: i64, e: i64) -> Vec<(i64, i64, bool, bool)> {
        let mut out = Vec::new();
        if self.end > self.n {
            let hi = e.min(self.end - self.n);
            if s < hi {
                out.push((s + self.n - self.start, hi + self.n - self.start, false, hi != e));
            }
        }
        let (lo, hi) = (s.max(self.start), e.min(self.end.min(self.n)));
        if lo < hi {
            out.push((lo - self.start, hi - self.start, lo != s, hi != e));
        }
        out
    }

    /// The flat list of pieces of `loc`, in order.
    fn parts(&self, loc: &Location) -> Vec<Location> {
        match loc {
            Location::Range((s, before), (e, after)) => self
                .pieces(*s, *e)
                .into_iter()
                // A "<" / ">" marker stays only on an end that was not cut.
                .map(|(a, b, cut_a, cut_b)| {
                    Location::Range((a, Before(before.0 && !cut_a)), (b, After(after.0 && !cut_b)))
                })
                .collect(),
            Location::Between(a, b) => {
                let inside = *a >= self.start && *a <= self.end.min(self.n);
                let after_origin = self.end > self.n && *a < self.end - self.n;
                if after_origin {
                    vec![Location::Between(a + self.n - self.start, b + self.n - self.start)]
                } else if inside {
                    vec![Location::Between(a - self.start, b - self.start)]
                } else {
                    Vec::new()
                }
            }
            Location::Complement(inner) => {
                let parts = self.parts(inner);
                if parts.is_empty() {
                    return parts;
                }
                // Same order inside as before: the complement reads the whole run backwards.
                vec![Location::Complement(Box::new(merged(parts, Location::Join)))]
            }
            Location::Join(v) | Location::Order(v) => v.iter().flat_map(|l| self.parts(l)).collect(),
            Location::Bond(v) => non_empty(v.iter().flat_map(|l| self.parts(l)).collect()).map(Location::Bond).into_iter().collect(),
            Location::OneOf(v) => non_empty(v.iter().flat_map(|l| self.parts(l)).collect()).map(Location::OneOf).into_iter().collect(),
            Location::External(..) | Location::Gap(_) => vec![loc.clone()],
        }
    }
}

/// `parts` as one location: touching neighbours merged, a single part as itself, else `wrap`.
fn merged(parts: Vec<Location>, wrap: fn(Vec<Location>) -> Location) -> Location {
    let mut out: Vec<Location> = Vec::new();
    for p in parts {
        // Merge only two plain ranges where the first ends exactly where the next starts.
        let touching = match (out.last(), &p) {
            (Some(Location::Range(_, (end, _))), Location::Range((s, _), _)) => end == s,
            _ => false,
        };
        if touching {
            if let (Some(Location::Range(from, _)), Location::Range(_, to)) = (out.pop(), p) {
                out.push(Location::Range(from, to));
            }
        } else {
            out.push(p);
        }
    }
    if out.len() == 1 { out.pop().unwrap() } else { wrap(out) }
}

/// Flatten a location into its base ranges, in order, plus the strand
/// (+1 forward, -1 reverse). Used for drawing and overlap queries.
pub fn spans(loc: &Location) -> (Vec<(i64, i64)>, i8) {
    let mut out = Vec::new();
    let mut strand = 1;
    collect(loc, false, &mut out, &mut strand);
    (out, strand)
}

// `&mut Vec<...>` lets the recursive calls append to one shared output list.
fn collect(loc: &Location, reverse: bool, out: &mut Vec<(i64, i64)>, strand: &mut i8) {
    match loc {
        Location::Range((s, _), (e, _)) => {
            out.push((*s, *e));
            if reverse {
                *strand = -1;
            }
        }
        // `20^21` is the gap between bases 20 and 21: at 0-based position 20, where the second base starts.
        Location::Between(_, b) => out.push((*b, *b)),
        Location::Complement(inner) => collect(inner, !reverse, out, strand),
        Location::Join(v) | Location::Order(v) | Location::Bond(v) | Location::OneOf(v) => {
            for l in v {
                collect(l, reverse, out, strand);
            }
        }
        Location::External(..) | Location::Gap(_) => {}
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn r(s: i64, e: i64) -> Location {
        Location::simple_range(s, e)
    }

    #[test]
    fn delete_rules() {
        assert_eq!(after_delete(&r(10, 20), 0, 5), Some(r(5, 15))); // before: shift
        assert_eq!(after_delete(&r(10, 20), 5, 15), Some(r(5, 10))); // overlaps start: trim
        assert_eq!(after_delete(&r(10, 20), 12, 14), Some(r(10, 18))); // inside: shrink
        assert_eq!(after_delete(&r(10, 20), 10, 20), None); // exactly covered: removed
        // complement(join(a,b)) losing one part collapses to complement(b)
        let cj = Location::Complement(Box::new(Location::Join(vec![r(0, 5), r(10, 20)])));
        assert_eq!(
            after_delete(&cj, 0, 5),
            Some(Location::Complement(Box::new(r(5, 15))))
        );
    }

    #[test]
    fn spans_and_strand() {
        let cj = Location::Complement(Box::new(Location::Join(vec![r(0, 5), r(10, 20)])));
        assert_eq!(spans(&cj), (vec![(0, 5), (10, 20)], -1));
        assert_eq!(spans(&r(3, 4)), (vec![(3, 4)], 1));
    }

    #[test]
    fn clip_window_keeps_parts_in_reading_order() {
        // join(191..200,1..10) on a 200-base circle, cut so the window starts at base 6:
        // old 190..199 -> 185..194 then old 0..4 -> 195..199 (touching: merged), then
        // old 5..9 -> 0..4. The old code put 0..4 before 195..199, scrambling the feature.
        let j = Location::Join(vec![r(190, 200), r(0, 10)]);
        let got = clip_window(&j, 5, 200, 200).unwrap();
        assert_eq!(spans(&got).0, vec![(185, 200), (0, 5)]);
        // A window that does not cut anything: just a shift.
        assert_eq!(clip_window(&r(10, 20), 5, 50, 200), Some(r(5, 15)));
        // Outside the window: gone. Partly inside: clipped.
        assert_eq!(clip_window(&r(100, 120), 5, 50, 200), None);
        assert_eq!(clip_window(&r(0, 20), 10, 50, 200), Some(r(0, 10)));
        // A complement keeps its strand and its parts' order.
        let c = Location::Complement(Box::new(Location::Join(vec![r(17, 20), r(0, 3)])));
        let got = clip_window(&c, 15, 10, 20).unwrap();
        assert_eq!(spans(&got), (vec![(2, 8)], -1));
        // A range split by the window start becomes two parts, the lower bases first.
        let got = clip_window(&r(8, 14), 10, 20, 20).unwrap();
        assert_eq!(spans(&got).0, vec![(18, 20), (0, 4)]);
        let got = clip_window(&r(8, 14), 12, 20, 20).unwrap();
        assert_eq!(spans(&got).0, vec![(16, 20), (0, 2)]);
    }
}
