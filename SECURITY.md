# Security policy

## Reporting a vulnerability

Please report security problems **privately**, not in a public issue.

1. Preferred: use GitHub's private vulnerability reporting. On the repository page open **Security > Report a vulnerability**
   (https://github.com/OpenDNA-project/openDNA/security/advisories/new). Only the maintainer sees it.
2. Or write to `<project email>`. Do not send sequence files you cannot share; describe how to reproduce instead.

Please include: the openDNA version (Help > About), what you did, what happened, and, if you can, a small file that triggers it.
This is a one-person, spare-time project: expect an acknowledgement within about a week and a fix or a clear answer as soon as it can be done. Please
give a reasonable time to fix a problem before you publish it; you will be credited in the release notes if you wish.

Examples of what is in scope: a crafted GenBank, FASTA, `.dna`, GFF / BED / VCF or BAM file that crashes the program or makes it run code or read
files it should not; the 3D structure page or the NCBI import doing something unsafe with data from a server; openDNA sending data you did not
expect (see [PRIVACY.md](PRIVACY.md)). Ordinary bugs (wrong result, crash with a normal file) are not vulnerabilities: use the issue templates.

## Supported versions

This is an early test version (0.1.x). Only the **latest release** is supported; fixes are not back-ported. Please update before reporting.

## What the project does and does not promise

- The program is offline by default and has no telemetry and no installer. It contacts the internet only for the features
  listed in [PRIVACY.md](PRIVACY.md).
- **Updates.** The exe can look for a newer release on GitHub and, when you press the button, replace itself with it. It only talks HTTPS to GitHub hosts,
  and refuses a download whose SHA-256 differs from the checksum file published with the release. That checksum comes from the same release page as the
  exe, so it catches a corrupted or truncated download or a tampering proxy, but **not a hijacked release** (an attacker would replace both files); only
  code signing fixes that, and the exe is not signed (next section). The check is off while the repository is private, and you can turn it off in Tools > Settings > Updates.
- File readers are written in Rust (GenBank via the `gb-io` crate, BAM / VCF via `noodles`, `.dna` with an XML reader) and are tested to give an
  error rather than crash on truncated and corrupt files, but they have **not been fuzzed or independently audited**. Treat files from untrusted sources
  with the care you would give any document.
- The *Local command* option of the 3D structure window runs whatever command you type, with your user rights. Do not paste commands you do not understand.
- The NCBI API key and email you may enter are stored in plain text in your data folder (see PRIVACY.md).
- Dependencies are not scanned automatically yet. The versions are in `core/Cargo.lock`, `requirements.txt` and [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## The Windows exe is not code-signed

The `oDNA_beta.exe` in the Releases is **unsigned** (a certificate costs money and is not in place for this test version), so Windows SmartScreen shows
"Windows protected your PC" the first time. That warning is expected. To make sure the file is the one that was published:

1. Download it only from https://github.com/OpenDNA-project/openDNA/releases (nowhere else; if you paid for it, you were scammed).
2. Each release lists a **SHA-256** checksum in its notes (and, when attached, in a `SHA256SUMS` file) for the file you download.
3. Compute the checksum of your download and compare it, letter for letter:

   ```powershell
   Get-FileHash .\oDNA_beta.exe -Algorithm SHA256
   ```

   or, in a command prompt, `certutil -hashfile oDNA_beta.exe SHA256`. If the release gives the checksum of a `.zip`, check the zip before you unpack it.
4. If the values differ, do not run the file, and report it as above.

The macOS and Linux downloads (`oDNA_beta-macos-arm64.zip`, `oDNA_beta-macos-x86_64.zip`, `oDNA_beta-linux-x86_64.tar.gz`) are not signed either: macOS refuses to
open the app the first time (right-click it, *Open*, then *Open* again, or System Settings > Privacy & Security > *Open Anyway*).
Check them the same way before you unpack them, with `shasum -a 256 <file>` (macOS) or `sha256sum <file>` (Linux), against the
`.sha256` file beside each download.

A matching checksum shows the file is the one that was published; it does not by itself prove the published file is harmless. The source is public, and you can build
the exe yourself (see the README).
