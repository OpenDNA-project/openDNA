# Rust primer for a Python user (with openDNA's own code)

Read alongside [core/src/fasta.rs](../core/src/fasta.rs) and [core/src/lib.rs](../core/src/lib.rs); every concept below appears there.

## 1. How the pieces fit

```
core/Cargo.toml          ← Rust project file (like pyproject.toml)
core/src/lib.rs          ← entry point; Python bindings (PyO3)
core/src/fasta.rs        ← plain Rust logic, no Python inside
core/pyproject.toml      ← tells pip to build the crate with maturin
app/main.py              ← Python GUI: `import opendna_core`
```

`maturin develop` compiles Rust into a native module (`opendna_core.pyd`) and installs it into your Python environment. Then Python imports it like any package. **Every time you change Rust code, rebuild** (`maturin develop -m core/Cargo.toml --release`). Python changes need no rebuild.

## 2. Commands you will use

| Command | Python equivalent | What it does |
|---|---|---|
| `cargo build` | – | compile (checks types, borrowing, everything) |
| `cargo test` | `pytest` | run the `#[test]` functions |
| `cargo fmt` | `black` | format code |
| `cargo clippy` | `ruff` | lint, suggests idiomatic fixes |
| `cargo add noodles-fasta` | `pip install` + edit pyproject | add a dependency to Cargo.toml |
| `maturin develop --release` | `pip install -e .` | build + install the Python module |

Always use `--release` when timing anything: debug builds can be 10–50× slower.

## 3. The language in ten ideas

**1. Static types, mostly inferred.** `let x = 5;` is an integer, and the compiler knows it. Function signatures are always explicit: `fn fetch(&mut self, name: &str, start: u64, end: u64) -> io::Result<Vec<u8>>`. Most bugs Python finds at runtime, Rust finds at compile time. A long compile error is normal; read the first one, it is usually precise.

**2. Variables are immutable unless `mut`.** `let mut pos = 0;` can change; `let e = ...;` cannot. Redeclaring (`let end = end.min(...)`) is allowed and called *shadowing*.

**3. `struct` + `impl` = class.** Data goes in `struct FaiEntry { ... }`, methods in `impl IndexedFasta { ... }`. No inheritance; shared behaviour uses *traits* (like Python protocols/ABCs). `#[derive(Debug, Clone)]` auto-implements common traits.

**4. Ownership: each value has exactly one owner.** When the owner goes out of scope, the value is freed. No garbage collector, no leaks, which is why Rust is fast and predictable with GB-sized data.

**5. Borrowing: `&` lends without giving away.**
- `&x` = read-only borrow; any number at once.
- `&mut x` = exclusive writable borrow; only one, and no readers meanwhile.

`fn fetch(&mut self, name: &str, ...)`: `self` is borrowed mutably (the file cursor moves), `name` is borrowed read-only. The compiler checks these rules, so data races and use-after-free can't happen. When the compiler complains about borrows, the usual fixes are: shorten how long a borrow lives, or `.clone()` the small value.

**6. `String` vs `&str`, `Vec<T>` vs `&[T]`.** `String`/`Vec` own their data (growable); `&str`/`&[T]` are borrowed views into someone else's data. Functions usually *take* the borrowed form and *return* the owned one.

**7. No `None` surprises: `Option<T>`.** A value that may be missing is `Some(v)` or `None`, and you must handle both: `match`, `if let Some(v) = ...`, `let Some(v) = ... else { ... }`, or `.unwrap_or(default)`.

**8. No exceptions: `Result<T, E>`.** A function that can fail returns `Ok(value)` or `Err(error)`. The `?` operator means "if Err, return it to my caller now", which is the equivalent of letting an exception propagate. At the Python boundary PyO3 turns an `Err` into a raised Python exception (`io::Error` becomes `OSError`). `.unwrap()` means "crash if Err"; fine in tests, avoid in app code.

**9. Closures and iterators** replace comprehensions:
```rust
self.inner.entries.iter().map(|e| (e.name.clone(), e.length)).collect()
// Python: [(e.name, e.length) for e in self.entries]
```
They compile to the same machine code as a hand-written loop, so there's no speed penalty.

**10. Expressions return values.** `if`/`match`/blocks produce values; the last line of a function **without `;`** is its return value.

## 4. Integers you will see

`u8` byte (one DNA base), `u32`/`u64` unsigned 32/64-bit, `i64` signed, `usize` the size type used for indexing. Positions are `u64` because genomes exceed 4 Gb. Converting between them is explicit: `n as usize`.

## 5. PyO3 in four macros

| Macro | Meaning |
|---|---|
| `#[pyclass]` | this struct becomes a Python class |
| `#[pymethods]` | methods in this `impl` are callable from Python |
| `#[new]` | the constructor (`__init__`) |
| `#[pymodule]` | the module function, run on `import` |

Types convert automatically: `String`/`&str` ↔ `str`, `u64` ↔ `int`, `Vec<T>` → `list`, tuples ↔ tuples, `PathBuf` ← `str`/`pathlib.Path`. `py.allow_threads(|| ...)` releases Python's GIL while Rust works, so the GUI keeps running.

**Design rule for this project:** Python never holds a whole genome. It asks Rust for a window (`fetch(name, start, end)`) and gets back only what's on screen.

## 6. Reading compiler errors

1. Read only the **first** error; later ones are often consequences.
2. The `help:` line usually contains the fix.
3. `rustc --explain E0502` prints a full explanation of any error code.
4. Borrow errors (E0499, E0502, E0505) mean two parts of the code want the same data at once; restructure or clone.

## 7. Further reading
- *The Rust Book*: https://doc.rust-lang.org/book/ (chapters 3–9 cover everything above)
- PyO3 guide: https://pyo3.rs
- Rust by Example: https://doc.rust-lang.org/rust-by-example/
