# Tutorial: Expose your local `ollama` LLM to the `CyteType` with `ngrok`

This tutorial will guide you through exposing your locally running `ollama` instance to the CyteType using `ngrok`. This method is incredibly fast and easy, making it perfect for temporary sharing, development, and testing.

- Pros: Extremely simple setup, no domain name required.
- Cons: The public URL is random and changes every time you restart `ngrok` (on the free plan). The connection dies when you close your terminal.

Prerequisites

- A Linux machine with ollama installed.
- At least one model pulled (e.g., `ollama pull qwen3:32b`).
- Access to the command line (CLI).

**NOTE**:

- Make sure the model you want to use supports *tools*. You can check out the [ollama library](https://ollama.com/library) for available models.
- Make sure you have enough resources available on your machine to run the model you want to use. You may encounter `Error 500: Internal Server Error` if the model is too large for your system's resources.

## Step 1: Configure Ollama to Listen on All Interfaces

By default, for security reasons, `ollama` only listens for requests from your local machine (`127.0.0.1`) at port `11434`. To allow a reverse proxy to connect to it, you need to configure it to listen on all network interfaces (`0.0.0.0`).

Create an override directory for the `ollama` service:

```bash
sudo mkdir -p /etc/systemd/system/ollama.service.d
```

Create a new configuration file:

```bash
sudo nano /etc/systemd/system/ollama.service.d/override.conf
```

Add the following content to the file:

```bash
[Service]
Environment="OLLAMA_HOST=0.0.0.0"
```

Save the file and exit the editor (in nano, press `Ctrl+O`, then press `Enter` to write, then `CTRL+X` to exit).

Reload the systemd daemon and restart `ollama`:

```bash
sudo systemctl daemon-reload
sudo systemctl restart ollama
```

Verify that Ollama is listening correctly:

Run this command and look for an entry showing `0.0.0.0:11434` or `*:11434`.

```bash
ss -tlnp | grep 11434
```

You should see output similar to:

```bash
LISTEN 0  4096       0.0.0.0:11434       0.0.0.0:*  users:(("ollama",pid=1234,fd=7))
```

Your `ollama` instance is now ready to accept connections from other devices on your local network and from our reverse proxy.

**NOTE:** If you are using WSL without systemd enabled you will have need to run the following on a *separate terminal*:

```bash
# OLLAMA_HOST=0.0.0.0 OLLAMA_KEEP_ALIVE=-1 ollama serve
OLLAMA_HOST=0.0.0.0 ollama serve
```

The `OLLAMA_KEEP_ALIVE` parameter controls how long the server will keep the model in memory after the response. By default, it's set to 5 minutes, but setting it to `-1` will keep the connection open indefinitely.

## Step 2: Set up and run `ngrok`

Install `ngrok`. The easiest way on Debian/Ubuntu is:

```bash
curl -sSL https://ngrok-agent.s3.amazonaws.com/ngrok.asc \
  | sudo tee /etc/apt/trusted.gpg.d/ngrok.asc >/dev/null \
  && echo "deb https://ngrok-agent.s3.amazonaws.com buster main" \
  | sudo tee /etc/apt/sources.list.d/ngrok.list \
  && sudo apt update \
  && sudo apt install ngrok
```

Or visit the `ngrok` [download page](https://ngrok.com/downloads/) for other installation methods.

## Step 3: Sign up and Authenticate

Go to the `ngrok` [dashboard](https://dashboard.ngrok.com/) and create a free account.

Find your authtoken on the dashboard.

Add the token to your local configuration (this is a one-time setup):

```bash
ngrok config add-authtoken <token>
```

## Step 4: Start the Tunnel & Get Your Public URL

This command tells `ngrok` to create a public URL and forward all traffic to your local port `11434`.

```bash
ngrok http 11434
```

`ngrok` will display a screen with your public URL. It will look something like this:

```bash
Session Status                online
Account                       Your Name (Plan: Free)
Version                       3.x.x
Region                        United States (us)
Forwarding                    https://random-string-of-letters.ngrok-free.app -> http://localhost:11434

Connections                   ttl     opn     rt1     rt5     p50     p90
                              0       0       0.00    0.00    0.00    0.00

```

**NOTE**: Keep this terminal window open! If you close it, the tunnel will stop.

## Step 5: Test Your Public OpenAI-Compatible API

You can now query your local model from anywhere using this public URL. The OpenAI-compatible endpoint is at `/v1`.

Example using `curl` (replace `https://your-ngrok-url.ngrok-free.app` with the URL from the previous step):

```bash
curl https://your-ngrok-url.ngrok-free.app/v1/chat/completions \
  -H "Content-Type: application/json" \
  -d '{
    "model": "qwen3:32b",
    "messages": [
      {
        "role": "user",
        "content": "Why is the sky blue?"
      }
    ]
  }'
```

Now you can use the following base url in your CyteType notebook as follows:

```python
# Run annotation
adata = annotator.run(
    study_context="...",
    model_config=[
        {
            "provider": "openai",
            "name": "qwen3:32b",
            "apiKey": "ollama",
            "baseUrl": "https://your-ngrok-url.ngrok-free.app/v1",
        }
    ]
)
```
