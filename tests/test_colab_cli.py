"""colab_cli: the credentials upload carries only the one netrc machine asked for, from a private temp file that does
not outlive the call, and no secret is ever part of code sent to the VM (the CLI logs executed code)."""
from pathlib import Path

from insar_wetlands import colab_cli as C


def test_netrc_upload_sends_one_machine_and_cleans_up(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    (home / ".netrc").write_text("machine other.example login a password SECRET1\n"
                                 "machine urs.earthdata.nasa.gov login u password SECRET2\n")
    monkeypatch.setattr(Path, "home", lambda: home)
    seen, sent_code = {}, []

    def fake_colab(*args, stdin=None, timeout=None):
        if args[0] == "upload":
            f = Path(args[3])
            seen["text"], seen["mode"], seen["path"] = f.read_text(), f.stat().st_mode & 0o777, f
        return 0

    monkeypatch.setattr(C, "colab", fake_colab)
    monkeypatch.setattr(C, "vm_python", lambda code, session, timeout=0: sent_code.append(code) or 0)
    assert C.cmd_netrc("s", "urs.earthdata.nasa.gov") == 0
    assert seen["text"] == "machine urs.earthdata.nasa.gov login u password SECRET2\n"
    assert seen["mode"] == 0o600 and not seen["path"].exists()
    assert all("SECRET" not in c for c in sent_code)
