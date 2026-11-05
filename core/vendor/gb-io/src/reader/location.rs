/// Parsers for location specifiers
use std::str;

use branch::alt;
use bytes::complete::{is_not, tag};
use character::complete::{char, i64};
use combinator::{map, map_opt, map_res, opt};
use multi::separated_list1;
use nom::combinator::{eof, peek, value};
use nom::{IResult, Parser};
use nom::{branch, bytes, character, combinator, multi};
use nom::sequence::{preceded, separated_pair, terminated};

use crate::seq::{After, Before, GapLength, Location};

/// openDNA patch. A location must end where a location can end: at a comma or a closing bracket
/// (inside join, complement...), at white space, or at the end of the text. Upstream stopped at the first
/// thing it did not understand and handed back what it had, so `10..(20.30)` (a fuzzy end) became the
/// single base `10` and the rest was silently dropped. Now such a text is not understood as a whole and
/// falls through to `location_external`, which keeps it word for word.
fn bounded<'a, P>(parser: P) -> impl Parser<&'a [u8], Output = Location, Error = nom::error::Error<&'a [u8]>>
where
    P: Parser<&'a [u8], Output = Location, Error = nom::error::Error<&'a [u8]>>,
{
    terminated(
        parser,
        peek(alt((tag(","), tag(")"), tag(" "), tag("\t"), tag("\r"), tag("\n"), map(eof, |_| &b""[..])))),
    )
}

fn location_single(input: &[u8]) -> IResult<&[u8], Location> {
    map(i64, |i| {
        Location::Range(((i - 1), Before(false)), (i, After(false)))
    }).parse(input)
}

fn location_between(input: &[u8]) -> IResult<&[u8], Location> {
    map_opt(
        separated_pair(i64, char('^'), i64),
        // This check is not 100% foolproof, because if b == 1, a must be
        // the last nt of a circular sequence, but since we might not yet
        // know the sequence length, we can't check this yet
        |(a, b)| {
            ((a - b).abs() == 1 || ((a == 1) ^ (b == 1))).then_some(Location::Between(a - 1, b - 1))
        },
    ).parse(input)
}

fn location_span(input: &[u8]) -> IResult<&[u8], Location> {
    let parser = (opt(char('<')), i64, tag(".."), opt(char('>')), i64);
    map(parser, |(before, a, _, after, b)| {
        Location::Range(
            (a - 1, Before(before.is_some())),
            (b, After(after.is_some())),
        )
    }).parse(input)
}

fn operator<'a, F, O>(name: &'static str, inner: F) -> impl Parser<&'a [u8], Output = O, Error = nom::error::Error<&'a [u8]>> 
where
    F: Parser<&'a [u8], Output = O, Error = nom::error::Error<&'a [u8]>>
{
    terminated(
        preceded(tag(name), preceded(char('('), inner)),
        char(')'),
    )
}

fn location_compound(input: &[u8]) -> IResult<&[u8], Location> {
    fn parser<'a>(
        name: &'static str,
        constructor: impl Fn(Vec<Location>) -> Location,
    ) -> impl Parser<&'a [u8], Output = Location, Error = nom::error::Error<&'a [u8]>> {
        map(
            operator(name, separated_list1(char(','), location)),
            constructor,
        )
    }
    alt((
        parser("join", Location::Join),
        parser("order", Location::Order),
        parser("bond", Location::Bond),
        parser("one-of", Location::OneOf),
    )).parse(input)
}

fn location_complement(input: &[u8]) -> IResult<&[u8], Location> {
    map(operator("complement", location), |l| {
        Location::Complement(Box::new(l))
    }).parse(input)
}

fn location_gap(input: &[u8]) -> IResult<&[u8], Location> {
    let gap_length = alt((
        map(i64, GapLength::Known),
        value(GapLength::Unk100, tag("unk100")),
        value(GapLength::Unknown, tag("")),
    ));
    map(operator("gap", gap_length), Location::Gap).parse(input)
}

fn location_external(input: &[u8]) -> IResult<&[u8], Location> {
    let accession_parser = map_res(is_not(": \t\r\n"), str::from_utf8);
    let location = opt(preceded(
        tag(":"),
        alt((
            // everything except location_external
            bounded(location_span),
            bounded(location_compound),
            bounded(location_complement),
            bounded(location_between),
            bounded(location_single),
            bounded(location_gap),
        )),
    ));
    map((accession_parser, location), |(accession, loc)| {
        Location::External(accession.to_owned(), loc.map(Box::new))
    }).parse(input)
}

pub fn location(input: &[u8]) -> IResult<&[u8], Location> {
    alt((
        bounded(location_span),
        bounded(location_compound),
        bounded(location_complement),
        bounded(location_between),
        bounded(location_single),
        bounded(location_gap),
        location_external,
    )).parse(input)
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_pos_single() {
        let s = "1\n";
        let p = location(s.as_bytes());
        match p {
            Ok((_, p)) => match p {
                Location::Range((0, Before(false)), (1, After(false))) => {}
                x => {
                    panic!("{:?}", x);
                }
            },
            x => {
                panic!("{:?}", x);
            }
        }
    }

    #[test]
    fn test_complement_join() {
        let p = location(&b"join(complement(2..6),8)"[..]).unwrap();
        match p {
            (_, Location::Join(inner)) => match &*inner {
                [Location::Complement(inner), Location::Range((7, Before(false)), (8, After(false)))] => {
                    match inner.as_ref() {
                        Location::Range((1, Before(false)), (6, After(false))) => {}
                        x => {
                            panic!("{:?}", x);
                        }
                    }
                }
                x => {
                    panic!("{:?}", x);
                }
            },
            x => {
                panic!("{:?}", x);
            }
        }
    }

    #[test]
    fn test_range() {
        let s = "44..>579";
        assert_eq!(
            location(s.as_bytes()),
            Ok((
                (&b""[..]),
                Location::Range((43, Before(false)), (579, After(true)))
            ))
        )
    }

    #[test]
    fn test_pos_join() {
        let s = "join(1,3)";
        let p = location(s.as_bytes());
        let inner = vec![Location::simple_range(0, 1), Location::simple_range(2, 3)];
        assert_eq!(p, Ok(((&b""[..]), Location::Join(inner))));
    }
}
