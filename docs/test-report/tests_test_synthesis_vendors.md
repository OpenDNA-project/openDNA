# `tests/test_synthesis_vendors.py`

[← all test files](README.md) · source: [`tests/test_synthesis_vendors.py`](../../tests/test_synthesis_vendors.py)

## `test_the_idt_login_and_screening_requests_have_the_shape_the_open_client_uses`

**Result:** ✅ passed in 0.00 s · [line 23](../../tests/test_synthesis_vendors.py#L23) · **Area:** Synthesis tab: company checks and stored logins · **Checked against:** hand-computed · **Strength:** 🔴 weak

**Why it matters:** A wrong login or screening request would send credentials to IDT in a form it rejects, or read the score from the wrong place.

Screening a sequence makes two calls: a password-grant login (client id and secret as HTTP Basic) and a screening call with the bearer token and the sequence as a JSON list; the score is the sum of the findings' scores.

| # | Input | Expected output |
|---|---|---|
| 1 | a 200-base sequence and fake credentials | token call to the identity server with grant_type=password, scope=test, the username and password, Basic base64(client id:secret); then a Bearer call whose body is [{Name, Sequence}]; score 12.5 from findings 4.5 and 8 |

**Review notes:**

- The shape being checked is the one the open-source SBOL-utilities client is reported to use, not one read from IDT's documentation (behind a login) or seen on the live service: the first real call by someone with an IDT API account is the real test.

## `test_an_idt_answer_is_read_whatever_its_nesting_and_a_missing_score_is_said`

**Result:** ✅ passed in 0.00 s · [line 39](../../tests/test_synthesis_vendors.py#L39) · **Area:** Synthesis tab: company checks and stored logins · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** An answer in an unexpected shape must not be shown as a good score.

Findings are collected from any nesting of lists and dictionaries; a boolean is not a score; an answer with no Score gives none; text that is not JSON is refused.

| # | Input | Expected output |
|---|---|---|
| 1 | {"Score": 3} | 3 |
| 2 | [[{Score 1}], [{Score 2.5}, {Name}]] | 3.5 |
| 3 | [{Score: true}] | no score |
| 4 | <html>maintenance</html> | error: not JSON |

**Review notes:**

- The shape being checked is the one the open-source SBOL-utilities client is reported to use, not one read from IDT's documentation (behind a login) or seen on the live service: the first real call by someone with an IDT API account is the real test.

## `test_a_length_idt_does_not_cover_is_refused_before_anything_is_sent`

**Result:** ✅ passed in 0.00 s · [line 49](../../tests/test_synthesis_vendors.py#L49) · **Area:** Synthesis tab: company checks and stored logins · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** IDT answers 0 outside its length range, which would read as an excellent score.

Sequences under 125 or over 3000 bases are refused with nothing sent; 125 and 3000 are accepted.

| # | Input | Expected output |
|---|---|---|
| 1 | 124 and 3001 bases | error naming the range and 'Nothing was sent', post never called |
| 2 | 125 and 3000 bases | sent |

**Review notes:**

- The range itself (125 to 3000) is as reported by the SBOL-utilities documentation, not read from IDT.

## `test_idt_errors_say_what_happened_without_repeating_the_login`

**Result:** ✅ passed in 0.00 s · [line 58](../../tests/test_synthesis_vendors.py#L58) · **Area:** Synthesis tab: company checks and stored logins · **Checked against:** hand-computed · **Strength:** 🟢 strong

**Why it matters:** An error message could leak the password or client secret into a screenshot, a log or a bug report.

A 401 from IDT gives its code and message and a hint about the four login values but none of their contents; an unreachable host and an answer with no token get their own messages.

| # | Input | Expected output |
|---|---|---|
| 1 | 401 Unauthorized with body invalid_client | message has '401 Unauthorized', 'invalid_client', 'client secret' and none of the four values |
| 2 | no route to host | Could not reach www.idtdna.com: no route |
| 3 | login answer without access_token | error: without an access token |

## `test_logins_go_to_the_systems_secure_store_or_stay_in_memory_never_in_a_file`

**Result:** ✅ passed in 0.00 s · [line 102](../../tests/test_synthesis_vendors.py#L102) · **Area:** Synthesis tab: company checks and stored logins · **Checked against:** hand-computed · **Strength:** 🟡 adequate

**Why it matters:** A password written to a file, or lost silently, or a remembered login that cannot be forgotten.

With a working system store the login is saved and loaded from it and forgotten from it; with no store (fail backend or no keyring) or one that refuses writes, it is held in memory only (save says so); forgetting empties both.

| # | Input | Expected output |
|---|---|---|
| 1 | a working fake keyring (Windows) | save True, stored JSON equals the login, load equals it, forget empties |
| 2 | fail backend / no keyring / write error | save False, load still returns it this run, forget empties |
| 3 | fail backend | backend() is None, secure() False |

**Review notes:**

- Uses a fake keyring: the real Windows, macOS and Linux stores are not exercised, and neither is the frozen build's backend discovery.
