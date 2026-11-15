import argparse
import json
import re
import subprocess
import wave
from pathlib import Path

import imageio_ffmpeg
from PIL import Image, ImageDraw, ImageFont


SCENES = [
    ("The comparison", None),
    ("Delivery age", "age_support.png"),
    ("Pregnancy-level results", "headline.png"),
    ("Gene ranking", "enrichment.png"),
    ("Where the RNA appears", "atlas.png"),
    ("What the marker check found", "composition.png"),
    ("What to measure next", None),
    ("Data and acknowledgments", None),
]
SIZE = (1920, 1080)
INK = "#182D3D"
TEAL = "#005E68"
SECONDARY = "#4D6271"


def font(size, bold=False):
    name = "segoeuib.ttf" if bold else "segoeui.ttf"
    return ImageFont.truetype(f"C:/Windows/Fonts/{name}", size)


def card(scene, figures, path):
    title, filename = scene
    frame = Image.new("RGB", SIZE, "white")
    draw = ImageDraw.Draw(frame)
    draw.rectangle((0, 0, 1920, 86), fill=INK)
    draw.text((52, 22), "OPENBIO DATATHON   /   TRACK 1", font=font(29, True), fill="white")
    draw.text((1866, 25), title, font=font(26), fill="#DAE9ED", anchor="ra")
    if filename:
        figure = Image.open(figures / filename).convert("RGB")
        figure.thumbnail((1820, 945), Image.Resampling.LANCZOS)
        frame.paste(figure, ((1920 - figure.width) // 2, 90 + (970 - figure.height) // 2))
    else:
        draw.rectangle((90, 235, 108, 794), fill=TEAL)
        if title == "The comparison":
            draw.text((152, 240), "Placental RNA when\npreeclampsia and growth\nrestriction occur together",
                      font=font(87, True), spacing=10, fill=INK)
            draw.text((155, 670), "A comparison of BMP-related gene expression at delivery",
                      font=font(43), fill=SECONDARY)
        elif title == "What to measure next":
            draw.text((152, 245), "Measure signaling inside cells", font=font(77, True), fill=INK)
            for row, text in enumerate([
                "Compare the two FGR groups at similar delivery ages.",
                "Measure endothelial BMPR2 RNA and nuclear pSMAD.",
                "Count the relative abundance of the same cell types.",
            ]):
                y = 410 + row * 112
                draw.ellipse((158, y + 17, 179, y + 38), fill=TEAL)
                draw.text((210, y - 3), text, font=font(43), fill=INK)
        else:
            draw.text((152, 240), "Data and acknowledgments", font=font(77, True), fill=INK)
            draw.text((154, 420), "Bulk RNA  •  OpenBio / Ackerman and Awamleh GSE114691\n"
                      "Cell reference  •  Admati and Skarbianskis\n"
                      "Pathway definitions  •  MSigDB Reactome 2024.1.Hs",
                      font=font(41), spacing=27, fill=INK)
            draw.text((154, 753), "Sources and rerun instructions are in the project README.",
                      font=font(38), fill=SECONDARY)
    frame.save(path)


def video(args):
    args.cache.mkdir(parents=True, exist_ok=True)
    prose = args.script.read_text(encoding="utf-8")
    narration = dict(re.findall(r"^### (.+?)\n\n(.+?)\s*(?=\n### |\Z)", prose, re.M | re.S))
    plan = {f"{i:02}": " ".join(narration[title].split()) for i, (title, _) in enumerate(SCENES, 1)}
    wav_plan = args.cache / "voice.json"
    wav_plan.write_text(json.dumps(plan), encoding="utf-8")
    ps1 = args.cache / "record_voice.ps1"
    ps1.write_text(r"""
param([string]$Plan, [string]$Folder)
Add-Type -AssemblyName System.Speech
$Voice = [System.Speech.Synthesis.SpeechSynthesizer]::new()
$Voice.SelectVoice('Microsoft Zira Desktop')
$Voice.Rate = -1
$Scenes = (Get-Content -LiteralPath $Plan -Raw -Encoding UTF8 | ConvertFrom-Json)
foreach ($Scene in $Scenes.PSObject.Properties) {
    $Voice.SetOutputToWaveFile((Join-Path $Folder "$($Scene.Name).wav"))
    $Voice.Speak([string] $Scene.Value)
    $Voice.SetOutputToNull()
}
$Voice.Dispose()
""", encoding="utf-8")
    subprocess.run(["powershell.exe", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File", str(ps1.resolve()),
                    str(wav_plan.resolve()), str(args.cache.resolve())], check=True)
    ffmpeg = imageio_ffmpeg.get_ffmpeg_exe()
    durations = []
    for i, scene in enumerate(SCENES, 1):
        slug = f"{i:02}"
        still = args.cache / f"{slug}.png"
        audio = args.cache / f"{slug}.wav"
        card(scene, args.figures, still)
        with wave.open(str(audio)) as sample:
            duration = sample.getnframes() / sample.getframerate() + .8
        durations.append(duration)
        segment = args.cache / f"{slug}.mp4"
        subprocess.run([
            ffmpeg, "-hide_banner", "-loglevel", "error", "-y", "-loop", "1", "-framerate", "24",
            "-i", str(still), "-i", str(audio), "-af", "apad", "-t", f"{duration:.3f}",
            "-c:v", "libx264", "-preset", "veryfast", "-tune", "stillimage", "-crf", "18",
            "-pix_fmt", "yuv420p", "-c:a", "aac", "-b:a", "128k", "-ar", "48000", str(segment),
        ], check=True)
    playlist = args.cache / "segments.txt"
    playlist.write_text("".join(f"file '{i:02}.mp4'\n" for i in range(1, len(SCENES) + 1)), encoding="ascii")
    if sum(durations) > 300:
        raise ValueError("Recording exceeds five minutes")
    subprocess.run([ffmpeg, "-hide_banner", "-loglevel", "error", "-y", "-f", "concat", "-safe", "0",
                    "-i", str(playlist), "-c", "copy", "-movflags", "+faststart", str(args.output)], check=True)
    elapsed = 0
    for (name, _), seconds in zip(SCENES, durations):
        print(f"{elapsed:.1f} s | {name}")
        elapsed += seconds
    print(f"{elapsed:.1f} s total")


if __name__ == "__main__":
    command = argparse.ArgumentParser()
    for argument in ("script", "figures", "cache", "output"):
        command.add_argument(f"--{argument}", required=True, type=Path)
    video(command.parse_args())
